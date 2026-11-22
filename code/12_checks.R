# ==============================================================================
# 12_checks.R  --  extra checks asked for in the internal review
#   A. residual variance: sigma^2-hat / sigma^2, and an alternative estimator
#      based on the variance of the ensemble's predictions under Q*
#   B. tail rule: pad of 0.10, 0.25 (default) and 0.50 x IQR of the split points
#   C. grid Monte Carlo error: 10 grids per cohort
#   (same cohorts and ensembles as Study 1, n = 1500, first 50 replicates)
#   D. PhysioNet: tail rule
# usage: Rscript 12_checks.R [workers]
# ==============================================================================

suppressPackageStartupMessages({ library(parallel); library(data.table) })
root_dir <- normalizePath(Sys.getenv("TREEIPF_ROOT", "."), winslash = "/")
sim_file <- file.path(root_dir, "code", "01_run_simulations.R")
source(sim_file)
args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args) >= 1) as.integer(args[1]) else 2L
out_dir <- file.path(root_dir, "results")

tasks_c <- tasks_for(1)[n == 1500L & rho %in% c(0, 0.7) & rep <= 50]

check_one <- function(i, tasks) {
  task <- as.list(tasks[i])
  set.seed(task$seed)
  X <- draw_x(task$n, task$rho, 0L); y <- draw_y(X, task$family, "linear", 0)
  gbt <- gbt_for(task, X, y, max(8L, round(task$n / 60)))
  rows <- list()
  for (pad in c(0.10, 0.25, 0.50)) {
    options(treeipf.pad = pad)
    set.seed(task$seed + 1L)
    f <- tree_ipf(gbt, M = 2200L, alpha = 0.05)
    s2_alt <- NA_real_
    if (task$family == "gaussian") {
      vm <- sum(f$q * (f$mu - sum(f$q * f$mu))^2)
      s2_alt <- (gbt$tss - gbt$n_hat * vm) / (gbt$n_hat - 6)
    }
    rows[[length(rows) + 1L]] <- data.table(check = "pad", pad = pad, grid = 0L, term = names(f$est),
                                            est = f$est, se = f$se, sigma2 = f$sigma2, sigma2_alt = s2_alt,
                                            rec_sd_x1 = f$moments$sd[1], rec_sd_x3 = f$moments$sd[3])
  }
  options(treeipf.pad = 0.25)
  for (g in 1:10) {
    set.seed(task$seed + 100L + g)
    f <- tree_ipf(gbt, M = 2200L, alpha = 0.05)
    rows[[length(rows) + 1L]] <- data.table(check = "grid", pad = 0.25, grid = g, term = names(f$est),
                                            est = f$est, se = f$se, sigma2 = f$sigma2, sigma2_alt = NA_real_,
                                            rec_sd_x1 = NA_real_, rec_sd_x3 = NA_real_)
  }
  d <- rbindlist(rows)
  d[, `:=`(family = task$family, rho = task$rho, rep = task$rep)]
  d
}

if (sys.nframe() == 0L) {
  cl <- makeCluster(workers)
  clusterExport(cl, c("sim_file"))
  clusterEvalQ(cl, { suppressPackageStartupMessages(library(data.table)); source(sim_file); NULL })
  clusterExport(cl, c("check_one", "tasks_c"))
  res <- rbindlist(parLapplyLB(cl, seq_len(nrow(tasks_c)), check_one, tasks = tasks_c))
  stopCluster(cl)
  fwrite(res, file.path(out_dir, "checks_sim_raw.csv"))

  # summaries
  pad_s <- res[check == "pad" & term != "(Intercept)", .(mean_est = mean(est), mean_se = mean(se)),
               by = .(family, rho, pad, term)]
  s2 <- res[check == "pad" & term == "X1" & family == "gaussian",
            .(sigma2_ratio = mean(sigma2) / 2.25, sigma2_alt_ratio = mean(sigma2_alt) / 2.25,
              rec_sd_x1 = mean(rec_sd_x1), rec_sd_x3 = mean(rec_sd_x3)), by = .(rho, pad)]
  icpt <- res[check == "pad" & term == "(Intercept)", .(mean_icpt = mean(est)), by = .(family, rho, pad)]
  gm <- res[check == "grid" & term != "(Intercept)", .(grid_sd = sd(est), se = mean(se)), by = .(family, rho, rep, term)]
  gm <- gm[, .(grid_sd = mean(grid_sd), mean_se = mean(se), ratio = mean(grid_sd / se)), by = .(family, rho, term)]
  fwrite(pad_s, file.path(out_dir, "checks_pad.csv")); fwrite(s2, file.path(out_dir, "checks_sigma2.csv"))
  fwrite(icpt, file.path(out_dir, "checks_intercept.csv")); fwrite(gm, file.path(out_dir, "checks_grid.csv"))

  # D. PhysioNet tail rule (same ensemble as the PhysioNet extras, seed 1001)
  suppressPackageStartupMessages(library(arrow))
  raw <- as.data.table(read_parquet(file.path(root_dir, "data", "physionet_ihd.parquet")))
  v <- c("age", "female", "hr_max", "sbp_min", "temp_max", "gcs_min", "bun_max", "creat_max", "wbc_max", "mechvent")
  s <- na.omit(raw[, c("y", v), with = FALSE])
  X <- as.matrix(data.frame(age10 = s$age / 10, female = s$female, hr10 = s$hr_max / 10, sbp10 = s$sbp_min / 10,
                            temp_max = s$temp_max, gcs_min = s$gcs_min, bun10 = s$bun_max / 10,
                            creat_max = s$creat_max, wbc10 = s$wbc_max / 10, mechvent = s$mechvent))
  set.seed(1001L)
  gbt <- fit_cv_gbt(X, as.integer(s$y), "binomial", B_max = 160L, lr = 0.06, maxdepth = 3L,
                    minbucket = 20L, subsample = 0.80, lambda = 2.0, nfolds = 5L)
  ph <- rbindlist(lapply(c(0.10, 0.25, 0.50), function(pad) {
    options(treeipf.pad = pad); set.seed(4343L)
    f <- tree_ipf(gbt, M = 5000L, alpha = 0.02)
    data.table(pad = pad, term = names(f$est), est = f$est, se = f$se,
               rec_sd = c(NA, f$moments$sd), emp_sd = c(NA, apply(X, 2, sd)))
  }))
  fwrite(ph, file.path(out_dir, "checks_pad_physionet.csv"))
  cat("checks done\n")
}
