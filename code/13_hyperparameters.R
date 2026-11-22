# ==============================================================================
# 13_hyperparameters.R  --  how much the tuning of the shared ensemble matters
#   part 1: same cohorts (Study 1, n = 1500, rho = 0.4), two other tunings
#           "fast": learning rate 0.3, 100 trees, depth 6, min leaf 5 (close to XGBoost defaults)
#           "slow": learning rate 0.02, up to 400 trees chosen by CV, depth 3
#   part 2: bootstrap when the analyst gets the hyperparameters wrong
#           (logistic, 50 cohorts, 20 refits): correct; subsampling assumed 1;
#           leaf penalty assumed 0; half the number of trees
# usage: Rscript 13_hyperparameters.R [workers]
# ==============================================================================

suppressPackageStartupMessages({ library(parallel); library(data.table) })
root_dir <- normalizePath(Sys.getenv("TREEIPF_ROOT", "."), winslash = "/")
sim_file <- file.path(root_dir, "code", "01_run_simulations.R")
source(sim_file)
args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args) >= 1) as.integer(args[1]) else 8L
rep_dir <- file.path(root_dir, "results", "reps", "S8")
dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)

base <- tasks_for(1)[n == 1500L & rho == 0.4]
jobs <- rbind(base[rep <= 100, .(family, seed, rep, part = "tuning")],
              base[family == "binomial" & rep <= 50, .(family, seed, rep, part = "wrong")])

job_one <- function(i, jobs) {
  jb <- as.list(jobs[i])
  f <- file.path(rep_dir, sprintf("%s_%s_%03d.rds", jb$part, jb$family, jb$rep))
  if (file.exists(f)) return(f)
  res <- tryCatch({
    set.seed(jb$seed)
    X <- draw_x(1500L, 0.4, 0L); y <- draw_y(X, jb$family, "linear", 0)
    target <- setNames(c(BETA[[jb$family]][1:5], 0), c("(Intercept)", colnames(X)))
    lam <- if (jb$family == "gaussian") 1 else 2
    out <- list()
    if (jb$part == "tuning") {
      fits <- list(
        fast = fit_cv_gbt(X, y, jb$family, lr = 0.3, maxdepth = 6L, minbucket = 5L, subsample = 0.8, lambda = lam, fixed_B = 100L),
        slow = fit_cv_gbt(X, y, jb$family, B_max = 400L, lr = 0.02, maxdepth = 3L, minbucket = 25L, subsample = 0.8, lambda = lam, nfolds = 3L))
      for (k in names(fits)) {
        set.seed(jb$seed + 7L)
        ft <- tree_ipf(fits[[k]], M = 2200L)
        out[[k]] <- data.table(arm = k, term = names(ft$est), est = ft$est, se = ft$se,
                               ens = fhat_on_true_x(fits[[k]], X, ft$cuts), B = length(fits[[k]]$trees))
      }
    } else {
      gbt <- fit_cv_gbt(X, y, "binomial", B_max = 110L, lr = 0.07, maxdepth = 3L, minbucket = 25L,
                        subsample = 0.8, lambda = 2, nfolds = 3L)
      arms <- list(correct = gbt, tau1 = gbt, lambda0 = gbt, halfB = gbt)
      arms$tau1$hyper$subsample <- 1; arms$tau1$n_hat <- round(gbt$trees[[1]]$size[1])
      arms$lambda0$hyper$lambda <- 0
      arms$halfB$hyper$B <- max(5L, round(gbt$hyper$B / 2))
      for (k in names(arms)) {
        set.seed(jb$seed + 7L)
        ft <- tree_ipf(arms[[k]], M = 2200L)
        set.seed(jb$seed + 11L)
        bt <- tree_ipf_boot(arms[[k]], ft, B_boot = 20L, M = 1500L)
        out[[k]] <- data.table(arm = k, term = names(ft$est), est = ft$est, se = ft$se,
                               se_boot = bt$se, bias_boot = bt$bias)
      }
    }
    d <- rbindlist(out, fill = TRUE)
    d[, `:=`(part = jb$part, family = jb$family, rep = jb$rep, target = target[term])]
    d
  }, error = function(e) list(error = conditionMessage(e)))
  saveRDS(res, f)
  f
}

if (sys.nframe() == 0L) {
  cl <- makeCluster(workers)
  clusterExport(cl, c("sim_file", "rep_dir"))
  clusterEvalQ(cl, { suppressPackageStartupMessages(library(data.table)); source(sim_file); NULL })
  clusterExport(cl, c("job_one", "jobs"))
  invisible(parLapplyLB(cl, seq_len(nrow(jobs)), job_one, jobs = jobs))
  stopCluster(cl)
  r <- lapply(list.files(rep_dir, full.names = TRUE), readRDS)
  r <- rbindlist(r[vapply(r, is.data.table, TRUE)], fill = TRUE)[term != "(Intercept)"]
  cv <- function(lo, hi, t) mean(lo <= t & t <= hi, na.rm = TRUE)
  s1 <- r[part == "tuning", .(rel_bias = mean(est - target) / target[1], ens_rel = mean(ens - target) / target[1],
                               coverage = cv(est - 1.96 * se, est + 1.96 * se, target), B = mean(B)),
          by = .(family, arm, term)]
  s2 <- r[part == "wrong", { bc <- est - bias_boot; s <- pmax(se, se_boot)
          .(rel_bias = mean(est - target) / target[1], rel_bias_bc = mean(bc - target) / target[1],
            cov_model = cv(est - 1.96 * se, est + 1.96 * se, target), cov_bc = cv(bc - 1.96 * s, bc + 1.96 * s, target)) },
          by = .(arm, term)]
  fwrite(s1, file.path(root_dir, "results", "sim_S8_tuning.csv"))
  fwrite(s2, file.path(root_dir, "results", "sim_S8_wrong_hyper.csv"))
  cat("hyperparameter study done\n")
}
