# ==============================================================================
# 09_privacy_leakage.R  --  how much of the cohort a shared ensemble gives away
#
# For Pima and PhysioNet we grow boosted ensembles (100 rounds, no CV) over a grid
# of tree depth and minimum leaf size, then ask what TreeIPF rebuilds from the file:
#   - cohort size n (exact or not)
#   - means and SDs of every predictor
#   - all pairwise correlations
# Three release options: leaf counts kept (default rpart / XGBoost covers / sklearn),
# counts stripped (only splits and leaf values; TreeIPF falls back to a uniform design).
#
# leakage scores (1 = rebuilt exactly, 0 = no better than knowing only the split ranges):
#   L_cor  = 1 - MAE(correlations) / MAE(all correlations set to 0)
#   L_mean = 1 - MAE(standardised means) / MAE(standardised means of the uniform design)
# ==============================================================================

root_dir <- normalizePath(Sys.getenv("TREEIPF_ROOT", "."), winslash = "/")
lib_file <- file.path(root_dir, "code", "00_tree_ipf_lib.R")
source(lib_file)
suppressPackageStartupMessages({ library(data.table); library(parallel); library(arrow) })
args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args) >= 1) as.integer(args[1]) else 8L
out_dir <- file.path(root_dir, "results"); rep_dir <- file.path(out_dir, "reps", "privacy")
dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)

load_cohort <- function(name) {
  if (name == "Pima") {
    d <- rbind(MASS::Pima.tr, MASS::Pima.te)
    return(list(X = as.matrix(d[, setdiff(names(d), "type")]), y = as.integer(d$type == "Yes")))
  }
  raw <- as.data.table(read_parquet(file.path(root_dir, "data", "physionet_ihd.parquet")))
  v <- c("age", "female", "hr_max", "sbp_min", "temp_max", "gcs_min",
         "bun_max", "creat_max", "wbc_max", "mechvent")
  s <- na.omit(raw[, c("y", v), with = FALSE])
  list(X = as.matrix(s[, ..v]), y = as.integer(s$y))
}

leak_scores <- function(X, Z, q) {
  mo <- weighted_moments(Z, q)
  ec <- cor(X); up <- upper.tri(ec)
  em <- colMeans(X); es <- apply(X, 2, sd)
  rc <- mo$cor[up]
  c(cor_r = cor(ec[up], rc), cor_mae = mean(abs(ec[up] - rc)),
    L_cor = 1 - mean(abs(ec[up] - rc)) / mean(abs(ec[up])),
    mean_err = mean(abs(mo$mean - em) / es),
    sd_ratio_err = mean(abs(mo$sd / es - 1)))
}

grid <- CJ(cohort = c("Pima", "PhysioNet"), depth = c(1L, 2L, 3L, 4L, 6L),
           minleaf = c(5L, 10L, 25L, 50L, 100L, 200L), seed = 1:3)

one <- function(i) {
  g <- as.list(grid[i])
  f <- file.path(rep_dir, sprintf("%s_d%d_m%d_s%d.rds", g$cohort, g$depth, g$minleaf, g$seed))
  if (file.exists(f)) return(readRDS(f))
  co <- load_cohort(g$cohort)
  if (g$cohort == "Pima" && g$minleaf > 100) { saveRDS(NULL, f); return(NULL) }
  set.seed(2000L + 10L * g$seed + g$depth)
  fo <- fit_cv_gbt(co$X, co$y, "binomial", lr = 0.1, maxdepth = g$depth, minbucket = g$minleaf,
                   subsample = 0.8, lambda = 2, fixed_B = 100L)
  ft <- tree_ipf(fo, M = 5000L, alpha = 0.05)
  un <- design_variant(fo, ft, "uniform")
  base_mean <- mean(abs(colMeans(un$Z) - colMeans(co$X)) / apply(co$X, 2, sd))
  s_counts <- leak_scores(co$X, ft$Z, ft$q)
  s_none <- leak_scores(co$X, un$Z, un$q)
  res <- rbind(
    data.table(release = "leaf counts kept", t(s_counts)),
    data.table(release = "counts stripped", t(s_none)))
  res[, L_mean := 1 - mean_err / base_mean]
  res[, `:=`(cohort = g$cohort, depth = g$depth, minleaf = g$minleaf, seed = g$seed,
             n = nrow(co$X), n_hat = c(fo$n_hat, NA), leaves = mean(sapply(fo$trees, function(t) sum(t$feat == 0L))))]
  saveRDS(res, f)
  res
}

if (sys.nframe() == 0L) {
  cl <- makeCluster(workers)
  clusterExport(cl, c("lib_file", "root_dir", "rep_dir", "grid", "load_cohort", "leak_scores", "one"))
  clusterEvalQ(cl, { source(lib_file); suppressPackageStartupMessages({ library(data.table); library(arrow) }); NULL })
  invisible(parLapplyLB(cl, seq_len(nrow(grid)), one))
  stopCluster(cl)
  res <- rbindlist(lapply(seq_len(nrow(grid)), one))
  fwrite(res, file.path(out_dir, "privacy_leakage.csv"))
  cat("privacy study done\n")
}
