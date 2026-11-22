# ==============================================================================
# 02_run_clinical_studies.R  --  real-data studies
#
#   A. single-cohort audits: Pima (n = 532, p = 7) and PhysioNet 2012 ICU
#      (n = 11,607, p = 10), each repeated over 50 seeds (CV split + ensemble)
#   B. PhysioNet extras: alpha = 0.05 sensitivity, temperature subgroup,
#      mechanical-ventilation diagnostics, and an XGBoost JSON model
#   C. International Stroke Trial, 12 countries: every country analysed both ways
#      (full-data GLM and TreeIPF from its boosted ensemble), parametric bootstrap,
#      then the meta-analyses (all-12 from ensembles = primary test; hybrids,
#      swapped split, classical-only, meta-regression on model type)
#
# usage: Rscript 02_run_clinical_studies.R [workers]
# ==============================================================================

root_dir <- normalizePath(Sys.getenv("TREEIPF_ROOT", "."), winslash = "/")
lib_file <- file.path(root_dir, "code", "00_tree_ipf_lib.R")
source(lib_file)
suppressPackageStartupMessages({
  library(arrow); library(data.table); library(metafor); library(parallel)
})
data_dir <- file.path(root_dir, "data")
out_dir  <- file.path(root_dir, "results")
rep_dir  <- file.path(out_dir, "reps", "clinical")
dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args) >= 1) as.integer(args[1]) else 8L
N_SEEDS <- 50L

# ------------------------------------------------------------------------------
# data
# ------------------------------------------------------------------------------

load_pima <- function() {
  d <- rbind(MASS::Pima.tr, MASS::Pima.te)
  list(X = d[, setdiff(names(d), "type")], y = as.integer(d$type == "Yes"))
}

load_physionet <- function() {
  raw <- as.data.table(read_parquet(file.path(data_dir, "physionet_ihd.parquet")))
  v <- c("age", "female", "hr_max", "sbp_min", "temp_max", "gcs_min",
         "bun_max", "creat_max", "wbc_max", "mechvent")
  s <- na.omit(raw[, c("y", v), with = FALSE])
  X <- data.frame(age10 = s$age / 10, female = s$female, hr10 = s$hr_max / 10,
                  sbp10 = s$sbp_min / 10, temp_max = s$temp_max, gcs_min = s$gcs_min,
                  bun10 = s$bun_max / 10, creat_max = s$creat_max, wbc10 = s$wbc_max / 10,
                  mechvent = s$mechvent)
  list(X = X, y = as.integer(s$y))
}

AUDIT <- list(
  Pima = list(B_max = 150L, lr = 0.05, minbucket = 10L, M = 5000L, alpha = 0.05),
  PhysioNet = list(B_max = 160L, lr = 0.06, minbucket = 20L, M = 5000L, alpha = 0.02)
)

# ------------------------------------------------------------------------------
# A. one audit run (one seed)
# ------------------------------------------------------------------------------

audit_once <- function(name, seed) {
  f <- file.path(rep_dir, sprintf("%s_%03d.rds", name, seed))
  if (file.exists(f)) return(readRDS(f))
  dat <- if (name == "Pima") load_pima() else load_physionet()
  cfg <- AUDIT[[name]]
  X <- as.matrix(dat$X); y <- dat$y
  set.seed(1000L * (name == "PhysioNet") + seed)
  gl <- ipd_glm(X, y, "binomial")
  gbt <- fit_cv_gbt(X, y, "binomial", B_max = cfg$B_max, lr = cfg$lr, maxdepth = 3L,
                    minbucket = cfg$minbucket, subsample = 0.80, lambda = 2.0, nfolds = 5L)
  ft <- tree_ipf(gbt, M = cfg$M, alpha = cfg$alpha, effects = TRUE)
  ab <- ablation_grid(gbt, ft, cfg$alpha)
  meth <- list(TreeIPF = ft, `Naive-Unif` = ab[["uniform|full|grid"]],
               `Unif-rootn` = ab[["uniform|trim|rootn"]], `Marg-IPF` = ab[["marginal|trim|rootn"]])
  rows <- rbindlist(lapply(names(meth), function(m) data.table(
    cohort = name, seed = seed, method = m, term = names(gl$est),
    est = as.numeric(meth[[m]]$est), se = as.numeric(meth[[m]]$se))))
  rows <- rbind(rows, data.table(cohort = name, seed = seed, method = "glm", term = names(gl$est),
                                 est = gl$est, se = gl$se),
                data.table(cohort = name, seed = seed, method = c("ALE", "PD"),
                           term = rep(names(gl$est), each = 2),
                           est = as.vector(rbind(ft$effects$ale, ft$effects$pd)),
                           se = rep(as.numeric(ft$se), each = 2)))
  emp_cor <- cor(X); up <- upper.tri(emp_cor)
  mom <- data.table(cohort = name, seed = seed, term = colnames(X),
                    var_type = ft$grid$var_type, emp_mean = colMeans(X), rec_mean = ft$moments$mean,
                    emp_sd = apply(X, 2, sd), rec_sd = ft$moments$sd,
                    n_splits = ft$n_splits, routed_mass = ft$routed_mass)
  cors <- data.table(cohort = name, seed = seed, emp = emp_cor[up], rec = ft$moments$cor[up],
                     pair = outer(colnames(X), colnames(X), paste, sep = ":")[up])
  res <- list(coef = rows, moments = mom, cors = cors,
              meta = data.table(cohort = name, seed = seed, n = nrow(X), n_hat = gbt$n_hat,
                                B = length(gbt$trees), kl = ft$kl, w_alpha = ft$w_alpha))
  saveRDS(res, f)
  res
}

# ------------------------------------------------------------------------------
# C. IST: one country
# ------------------------------------------------------------------------------

IST_VARS <- c("age10", "female", "sbp10", "consc_drowsy", "ratrial", "ct_inf", "rdef1", "rdef2")

load_ist <- function() {
  ist <- as.data.table(read_parquet(file.path(data_dir, "ist_dd6.parquet")))
  map <- fread(file.path(data_dir, "ist_country_codes.csv"))
  ist <- merge(ist, map, by = "country", all.x = TRUE)
  ist[, `:=`(age10 = age / 10, sbp10 = rsbp / 10,
             consc_drowsy = as.integer(rconsc != "F"), ct_inf = as.integer(ct == "infarct"))]
  cc <- na.omit(ist[, c("country", "country_code", "y", IST_VARS), with = FALSE])
  top <- cc[, .(N = .N, events = sum(y), prev = mean(y)), by = .(country, country_code)][order(-N)][1:12]
  top[, rank_id := .I]
  list(cc = cc[country %in% top$country], top = top)
}

ist_country <- function(k, B_boot = 30L) {
  f <- file.path(rep_dir, sprintf("IST_%02d.rds", k))
  if (file.exists(f)) return(readRDS(f))
  L <- load_ist(); top <- L$top
  sub <- L$cc[country == top$country[k]]
  X <- as.matrix(sub[, ..IST_VARS]); y <- as.integer(sub$y)
  gl <- ipd_glm(X, y, "binomial")
  set.seed(303600L + k)
  gbt <- fit_cv_gbt(X, y, "binomial", B_max = 150L, lr = 0.06, maxdepth = 3L,
                    minbucket = max(10L, round(nrow(X) / 50)), subsample = 0.80,
                    lambda = 2.0, nfolds = 5L)
  ft <- tree_ipf(gbt, M = 4000L, alpha = 0.05)
  ab <- ablation_grid(gbt, ft, 0.05)
  bt <- tree_ipf_boot(gbt, ft, B_boot = B_boot, M = 2000L, alpha = 0.05)
  res <- data.table(rank_id = k, country_code = top$country_code[k], n = nrow(X),
                    n_hat = gbt$n_hat, B = length(gbt$trees), term = names(gl$est),
                    glm_est = gl$est, glm_se = gl$se,
                    ipf_est = ft$est, ipf_se = ft$se,
                    boot_se = bt$se, boot_bias = bt$bias,
                    unif_est = ab[["uniform|full|grid"]]$est, unif_se = ab[["uniform|full|grid"]]$se,
                    unifr_est = ab[["uniform|trim|rootn"]]$est, unifr_se = ab[["uniform|trim|rootn"]]$se,
                    marg_est = ab[["marginal|trim|rootn"]]$est, marg_se = ab[["marginal|trim|rootn"]]$se)
  saveRDS(res, f)
  res
}

# ------------------------------------------------------------------------------
# meta-analyses
# ------------------------------------------------------------------------------

pool <- function(yi, sei, label, tm, ml = NULL) {
  # REML; step halving if Fisher scoring struggles; DerSimonian-Laird only as a last resort
  m <- tryCatch(rma(yi = yi, sei = sei, method = "REML"), error = function(e)
    tryCatch(rma(yi = yi, sei = sei, method = "REML", control = list(stepadj = 0.5, maxiter = 1000)),
             error = function(e) rma(yi = yi, sei = sei, method = "DL")))
  w <- weights(m)
  data.table(term = tm, strategy = label, K = length(yi), pooled = as.numeric(m$beta),
             se = m$se, lb = m$ci.lb, ub = m$ci.ub, OR = exp(as.numeric(m$beta)),
             OR_lb = exp(m$ci.lb), OR_ub = exp(m$ci.ub), tau2 = m$tau2, I2 = m$I2,
             ml_weight = if (is.null(ml)) NA_real_ else sum(w[ml]), tau_method = m$method)
}

meta_all <- function(ce) {
  out <- list()
  for (tm in IST_VARS) {
    d <- ce[term == tm][order(rank_id)]
    odd <- d$rank_id %% 2L == 1L
    bc <- d$ipf_est - d$boot_bias                 # bootstrap bias-corrected estimate
    bse <- pmax(d$boot_se, d$ipf_se)
    out[[tm]] <- rbind(
      pool(d$glm_est, d$glm_se, "Gold: GLM in all 12", tm),
      pool(d$ipf_est, d$ipf_se, "All 12 from ensembles: TreeIPF", tm),
      pool(bc, bse, "All 12 from ensembles: TreeIPF + bootstrap", tm),
      pool(d$unif_est, d$unif_se, "All 12 from ensembles: Naive-Unif", tm),
      pool(d$unifr_est, d$unifr_se, "All 12 from ensembles: Unif-rootn", tm),
      pool(d$marg_est, d$marg_se, "All 12 from ensembles: Marg-IPF", tm),
      pool(ifelse(odd, d$ipf_est, d$glm_est), ifelse(odd, d$ipf_se, d$glm_se), "Hybrid (odd ML): TreeIPF", tm, odd),
      pool(ifelse(odd, d$unif_est, d$glm_est), ifelse(odd, d$unif_se, d$glm_se), "Hybrid (odd ML): Naive-Unif", tm, odd),
      pool(ifelse(!odd, d$ipf_est, d$glm_est), ifelse(!odd, d$ipf_se, d$glm_se), "Hybrid (even ML): TreeIPF", tm, !odd),
      pool(ifelse(!odd, d$unif_est, d$glm_est), ifelse(!odd, d$unif_se, d$glm_se), "Hybrid (even ML): Naive-Unif", tm, !odd),
      pool(d$glm_est[!odd], d$glm_se[!odd], "Classical only (6 even)", tm),
      pool(d$glm_est[odd], d$glm_se[odd], "Classical only (6 odd)", tm))
  }
  rbindlist(out)
}

meta_regression <- function(ce) {
  rbindlist(lapply(IST_VARS, function(tm) {
    d <- ce[term == tm][order(rank_id)]
    odd <- as.integer(d$rank_id %% 2L == 1L)
    m <- rma(yi = ifelse(odd == 1, d$ipf_est, d$glm_est),
             sei = ifelse(odd == 1, d$ipf_se, d$glm_se), mods = ~ odd, method = "REML")
    data.table(term = tm, ml_shift = m$beta[2], ml_shift_se = m$se[2], p = m$pval[2],
               I2_without = m$I2, QM_p = m$QMp)
  }))
}

# ------------------------------------------------------------------------------
# run
# ------------------------------------------------------------------------------

if (sys.nframe() == 0L) {
  cl <- makeCluster(workers)
  clusterExport(cl, c("lib_file", "root_dir", "data_dir", "rep_dir", "AUDIT", "IST_VARS",
                      "load_pima", "load_physionet", "audit_once", "load_ist", "ist_country"))
  clusterEvalQ(cl, { source(lib_file); suppressPackageStartupMessages({ library(data.table); library(arrow) }); NULL })
  jobs <- c(lapply(seq_len(12), function(k) list("IST", k)),
            lapply(seq_len(N_SEEDS), function(s) list("PhysioNet", s)),
            lapply(seq_len(N_SEEDS), function(s) list("Pima", s)))
  invisible(parLapplyLB(cl, jobs, function(j) {
    if (j[[1]] == "IST") ist_country(j[[2]]) else audit_once(j[[1]], j[[2]]); NULL
  }))
  stopCluster(cl)

  # ---- A: collect the audits
  au <- lapply(c("Pima", "PhysioNet"), function(nm) lapply(seq_len(N_SEEDS), function(s) audit_once(nm, s)))
  au <- unlist(au, recursive = FALSE)
  fwrite(rbindlist(lapply(au, `[[`, "coef")), file.path(out_dir, "clinical_audit_coefficients.csv"))
  fwrite(rbindlist(lapply(au, `[[`, "moments")), file.path(out_dir, "clinical_audit_moments.csv"))
  fwrite(rbindlist(lapply(au, `[[`, "cors")), file.path(out_dir, "clinical_audit_correlations.csv"))
  fwrite(rbindlist(lapply(au, `[[`, "meta")), file.path(out_dir, "clinical_audit_runs.csv"))

  # ---- B: PhysioNet extras (seed 1)
  phy <- load_physionet(); X <- as.matrix(phy$X); y <- phy$y
  set.seed(1001L)
  gbt <- fit_cv_gbt(X, y, "binomial", B_max = 160L, lr = 0.06, maxdepth = 3L,
                    minbucket = 20L, subsample = 0.80, lambda = 2.0, nfolds = 5L)
  f02 <- tree_ipf(gbt, M = 5000L, alpha = 0.02); f05 <- tree_ipf(gbt, M = 5000L, alpha = 0.05)
  gl <- ipd_glm(X, y, "binomial")
  unif <- design_variant(gbt, f02, "uniform")
  tq <- data.table(
    item = c("share temp_max <= 36.6", "mortality temp_max <= 36.6", "mortality 36.6 < temp_max <= 38.0",
             "uniform design share below 37.0", "TreeIPF share below 37.0", "data share below 37.0",
             "temp_max mean (data)", "temp_max sd (data)", "temp_max mean (TreeIPF)", "temp_max sd (TreeIPF)",
             "uniform design lower edge", "uniform design upper edge",
             "mechvent splits", "mechvent routed mass", "gcs_min splits", "cor(mechvent, gcs_min)"),
    value = c(mean(X[, "temp_max"] <= 36.6), mean(y[X[, "temp_max"] <= 36.6]),
              mean(y[X[, "temp_max"] > 36.6 & X[, "temp_max"] <= 38]),
              mean(unif$Z[, "temp_max"] < 37), sum(f02$q[f02$Z[, "temp_max"] < 37]), mean(X[, "temp_max"] < 37),
              mean(X[, "temp_max"]), sd(X[, "temp_max"]),
              f02$moments$mean["temp_max"], f02$moments$sd["temp_max"],
              min(unif$Z[, "temp_max"]), max(unif$Z[, "temp_max"]),
              f02$n_splits[10], f02$routed_mass[10], f02$n_splits[6], cor(X[, "mechvent"], X[, "gcs_min"])))
  fwrite(tq, file.path(out_dir, "physionet_details.csv"))
  sens <- data.table(term = names(gl$est), glm = gl$est, glm_se = gl$se,
                     ipf_a02 = f02$est, se_a02 = f02$se, ipf_a05 = f05$est, se_a05 = f05$se)
  fwrite(sens, file.path(out_dir, "physionet_alpha_sensitivity.csv"))

  # XGBoost model file on PhysioNet
  suppressPackageStartupMessages(library(xgboost))
  set.seed(4242L)
  par <- list(objective = "binary:logistic", eta = 0.1, max_depth = 3, subsample = 0.8,
              lambda = 1, min_child_weight = 2, nthread = 4)
  dm <- xgb.DMatrix(X, label = y)
  cv <- xgb.cv(par, dm, nrounds = 500, nfold = 5, early_stopping_rounds = 20, verbose = 0)
  best <- cv$early_stop$best_iteration; if (is.null(best)) best <- cv$best_iteration
  bst <- xgb.train(par, dm, nrounds = best, verbose = 0)
  xp <- file.path(out_dir, "libmodels", "physionet_xgboost.json")
  dir.create(dirname(xp), showWarnings = FALSE, recursive = TRUE)
  xgb.save(bst, xp)
  fx <- read_xgb_json(xp, subsample = 0.8)
  set.seed(4243L); tx <- tree_ipf(fx, M = 5000L, alpha = 0.02)
  set.seed(4243L); tc <- tree_ipf(fx, M = 5000L, alpha = 0.02, cover_mode = "count")
  fwrite(data.table(term = names(gl$est), glm = gl$est, glm_se = gl$se,
                    xgb_ipf = tx$est, xgb_se = tx$se, xgb_count = tc$est, xgb_count_se = tc$se,
                    n_hat = fx$n_hat, rounds = best),
         file.path(out_dir, "physionet_xgboost.csv"))

  # ---- C: IST meta-analyses
  L <- load_ist()
  fwrite(L$top, file.path(out_dir, "ist_top12_countries.csv"))
  ce <- rbindlist(lapply(seq_len(12), ist_country))
  fwrite(ce, file.path(out_dir, "ist_country_estimates.csv"))
  fwrite(meta_all(ce), file.path(out_dir, "ist_meta_analysis.csv"))
  fwrite(meta_regression(ce), file.path(out_dir, "ist_meta_regression.csv"))
  cat("clinical studies done\n")
}
