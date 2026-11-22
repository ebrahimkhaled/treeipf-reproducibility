# does running SMART longer than 150 updates change the estimates?
# 20 Study 1 cohorts per outcome type (n = 1500, rho = 0.4), same grid, 150 vs 400 updates
suppressPackageStartupMessages(library(data.table))
root_dir <- normalizePath(Sys.getenv("TREEIPF_ROOT", "."), winslash = "/")
source(file.path(root_dir, "code", "01_run_simulations.R"))
tk <- tasks_for(1)[n == 1500L & rho == 0.4 & rep <= 20]
res <- rbindlist(lapply(seq_len(nrow(tk)), function(i) {
  task <- as.list(tk[i]); set.seed(task$seed)
  X <- draw_x(1500L, 0.4, 0L); y <- draw_y(X, task$family, "linear", 0)
  g <- gbt_for(task, X, y, 25L)
  set.seed(9); a <- tree_ipf(g, max_iter = 150L)
  set.seed(9); b <- tree_ipf(g, max_iter = 400L)
  data.table(family = task$family, term = names(a$est)[-1], d = (b$est - a$est)[-1], se = a$se[-1], kl150 = a$kl, kl400 = b$kl)
}))
out <- res[, .(mean_abs_change = mean(abs(d), na.rm = TRUE), max_change_in_se = max(abs(d) / se, na.rm = TRUE),
               kl150 = mean(kl150), kl400 = mean(kl400)), by = family]
print(out)
fwrite(out, file.path(root_dir, "results", "checks_iterations.csv"))
