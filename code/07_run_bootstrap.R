# ==============================================================================
# 07_run_bootstrap.R  --  parametric bootstrap for TreeIPF (Study 6)
#
# For each simulated cohort we rebuild the same ensemble as in Study 1 / 5
# (same seed, same code path), run TreeIPF, then draw B_BOOT cohorts from
# (Q*, f-hat, sigma^2-hat), retrain with the released hyperparameters and rerun
# TreeIPF.  The spread gives a bootstrap SE; the mean shift estimates the bias
# that comes from retraining (mostly regularisation shrinkage).
#
# usage: Rscript 07_run_bootstrap.R [workers] [B_BOOT] [R]
# ==============================================================================

suppressPackageStartupMessages({ library(parallel); library(data.table) })
root_dir <- normalizePath(Sys.getenv("TREEIPF_ROOT", "."), winslash = "/")
sim_file <- file.path(root_dir, "code", "01_run_simulations.R")
source(sim_file)          # brings the library, the DGPs and tasks_for()

args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args) >= 1) as.integer(args[1]) else 8L
B_BOOT  <- if (length(args) >= 2) as.integer(args[2]) else 30L
R_BOOT  <- if (length(args) >= 3) as.integer(args[3]) else 100L

t1 <- tasks_for(1)[n %in% c(600L, 1500L) & rho %in% c(0, 0.7) & rep <= R_BOOT]
t5 <- tasks_for(5)[scenario == "nonlinear" & n == 1500L & rep <= R_BOOT]
tasks_b <- rbind(t1, t5)

boot_file <- function(task) {
  file.path(root_dir, "results", "reps", "S6",
            sprintf("%s_%s_n%d_r%.1f_%03d.rds", task$family, task$scenario, task$n, task$rho, task$rep))
}

boot_one <- function(i, tasks) {
  task <- as.list(tasks[i])
  f <- boot_file(task)
  if (file.exists(f)) return(f)
  res <- tryCatch({
    set.seed(task$seed)
    p_extra <- 0L
    X <- draw_x(task$n, task$rho, p_extra)
    y <- draw_y(X, task$family, task$scenario, task$b5)
    minb <- if (task$study == "S1") max(8L, round(task$n / 60)) else 10L
    gbt <- gbt_for(task, X, y, minb)
    fit <- tree_ipf(gbt, M = 2200L, alpha = task$alpha)
    theta_I <- if (task$scenario == "linear") {
      setNames(c(BETA[[task$family]][1:5], 0), names(fit$est))
    } else true_projection(task$family, task$scenario, task$rho, 0L, fit$cuts)
    set.seed(task$seed + 5000000L)
    bt <- tree_ipf_boot(gbt, fit, B_boot = B_BOOT, M = 1500L, alpha = task$alpha)
    data.table(study = "S6", family = task$family, scenario = task$scenario, n = task$n,
               rho = task$rho, rep = task$rep, term = names(fit$est),
               est = fit$est, se_model = fit$se, se_boot = bt$se, bias_boot = bt$bias,
               q025 = apply(bt$draws, 2, quantile, 0.025, na.rm = TRUE),
               q975 = apply(bt$draws, 2, quantile, 0.975, na.rm = TRUE),
               n_ok = sum(complete.cases(bt$draws[, -1, drop = FALSE])),
               target = as.numeric(theta_I[names(fit$est)]))
  }, error = function(e) list(error = conditionMessage(e)))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  saveRDS(res, f)
  f
}

if (sys.nframe() == 0L) {
  todo <- which(!file.exists(vapply(seq_len(nrow(tasks_b)), function(i) boot_file(as.list(tasks_b[i])), "")))
  cat(sprintf("[%s] bootstrap: %d cohorts, %d to run, B = %d, %d workers\n",
              format(Sys.time(), "%H:%M"), nrow(tasks_b), length(todo), B_BOOT, workers))
  if (length(todo)) {
    cl <- makeCluster(workers)
    clusterExport(cl, c("sim_file", "root_dir", "B_BOOT"))
    clusterEvalQ(cl, { suppressPackageStartupMessages(library(data.table)); source(sim_file); NULL })
    clusterExport(cl, c("boot_file", "boot_one", "tasks_b"))
    invisible(parLapplyLB(cl, todo, boot_one, tasks = tasks_b))
    stopCluster(cl)
  }
  cat(sprintf("[%s] bootstrap done\n", format(Sys.time(), "%H:%M")))
}
