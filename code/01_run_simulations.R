# ==============================================================================
# 01_run_simulations.R  --  simulation studies for the TreeIPF paper (ADEMP)
#
#   Study 1  main grid: n {300,600,1500} x rho {0,.4,.7} x family     (R = 200)
#            methods: oracles, TreeIPF (GBT, bagged), the 3x2x2 ablation,
#            ALE and PD average slopes
#   Study 2  tree depth d {1,2,3,5}, n = 800, rho = .7                 (R = 200)
#   Study 3  trimming budget alpha {0,.02,.05,.10,.20}, n = 800, rho = .4
#   Study 4  power: weak X5 effect, n = 600, rho = .4                   (R = 200)
#   Study 5  scope: nonlinear + interaction + heteroscedastic truth,
#            rare events (about 10%), p = 12                            (R = 200)
#
# usage:  Rscript 01_run_simulations.R [studies] [workers]
#         e.g. Rscript 01_run_simulations.R 1,2,3 8
# every replicate is written to results/reps/<study>/ and skipped on rerun
# ==============================================================================

suppressPackageStartupMessages({ library(parallel); library(data.table) })
root_dir <- normalizePath(Sys.getenv("TREEIPF_ROOT", "."), winslash = "/")
lib_file <- file.path(root_dir, "code", "00_tree_ipf_lib.R")
source(lib_file)

args <- commandArgs(trailingOnly = TRUE)
studies <- if (length(args) >= 1) as.integer(strsplit(args[1], ",")[[1]]) else 1:5
workers <- if (length(args) >= 2) as.integer(args[2]) else 8L
R_REP <- 200L

# ------------------------------------------------------------------------------
# data-generating mechanisms
# ------------------------------------------------------------------------------

BETA <- list(gaussian = c(1.0, 1.2, -0.9, 0.7, 0.8, 0.0),
             binomial = c(-0.4, 0.9, -0.7, 0.6, 0.8, 0.0))

draw_x <- function(n, rho, p_extra = 0L) {
  X12 <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, rho, rho, 1), 2))
  X <- cbind(X1 = X12[, 1], X2 = X12[, 2], X3 = rexp(n) - 1,
             X4 = rbinom(n, 1, 0.4), X5 = rnorm(n))
  if (p_extra > 0) {
    # X6, X7 correlated with X1 (weak effects), X8 binary, X9-X12 noise
    e <- cbind(X6 = 0.3 * X[, 1] + sqrt(1 - 0.09) * rnorm(n),
               X7 = 0.3 * X[, 1] + sqrt(1 - 0.09) * rnorm(n),
               X8 = rbinom(n, 1, 0.3),
               X9 = rnorm(n), X10 = rexp(n) - 1, X11 = rnorm(n), X12 = rbinom(n, 1, 0.5))
    X <- cbind(X, e)
  }
  X
}

# true conditional mean (link scale) for each scenario
true_eta <- function(X, family, scenario, b5 = 0) {
  b <- BETA[[family]]
  if (scenario == "linear") {
    return(b[1] + X[, 1:5] %*% c(b[2:5], b5))
  }
  if (scenario == "rare") {
    return(-2.9 + X[, 1:5] %*% b[2:6])
  }
  if (scenario == "p12") {
    extra <- if (family == "gaussian") c(0.3, -0.3, 0.4, 0, 0, 0, 0) else c(0.25, -0.25, 0.4, 0, 0, 0, 0)
    return(b[1] + X[, 1:5] %*% b[2:6] + X[, 6:12] %*% extra)
  }
  # nonlinear: J-shaped skewed biomarker and an X1 x X4 interaction
  x3 <- X[, 3]
  if (family == "gaussian") {
    1.0 + 1.2 * X[, 1] - 0.9 * X[, 2] + 0.45 * x3 + 0.30 * x3^2 + 0.8 * X[, 4] + 0.6 * X[, 1] * X[, 4]
  } else {
    -0.6 + 0.9 * X[, 1] - 0.7 * X[, 2] + 0.35 * x3 + 0.25 * x3^2 + 0.8 * X[, 4] + 0.5 * X[, 1] * X[, 4]
  }
}

draw_y <- function(X, family, scenario, b5 = 0) {
  eta <- as.vector(true_eta(X, family, scenario, b5))
  if (family == "binomial") return(rbinom(nrow(X), 1, sig(eta)))
  sdv <- if (scenario == "nonlinear") 1.5 * exp(0.25 * X[, 1]) / exp(0.25^2 / 2) else 1.5
  eta + rnorm(nrow(X), 0, sdv)
}

# working-model projections of the true mean over a large sample from P_X
BIG <- new.env()
big_sample <- function(rho, p_extra) {
  key <- paste(rho, p_extra)
  if (is.null(BIG[[key]])) {
    old <- .Random.seed; set.seed(99991)
    BIG[[key]] <- draw_x(200000L, rho, p_extra)
    .Random.seed <<- old
  }
  BIG[[key]]
}

true_projection <- function(family, scenario, rho, p_extra, cuts = NULL, b5 = 0) {
  Xb <- big_sample(rho, p_extra)
  rows <- if (is.null(cuts)) rep(TRUE, nrow(Xb)) else in_box(Xb, cuts)
  mu <- as.vector(true_eta(Xb[rows, , drop = FALSE], family, scenario, b5))
  if (family == "binomial") mu <- sig(mu)
  D <- cbind(1, Xb[rows, , drop = FALSE])
  setNames(wls_glm(D, rep(1 / sum(rows), sum(rows)), mu, family), c("(Intercept)", colnames(Xb)))
}

# ------------------------------------------------------------------------------
# one replicate
# ------------------------------------------------------------------------------

pack <- function(task, method, est, se, target, extra = list()) {
  d <- data.table(study = task$study, family = task$family, scenario = task$scenario,
                  n = task$n, rho = task$rho, depth = task$depth, alpha = task$alpha,
                  b5 = task$b5, rep = task$rep, method = method,
                  term = names(est), est = as.numeric(est), se = as.numeric(se),
                  target = as.numeric(target[names(est)]))
  for (k in names(extra)) d[[k]] <- extra[[k]]
  d
}

gbt_for <- function(task, X, y, minb) {
  fit_cv_gbt(X, y, task$family, B_max = 110L, lr = 0.07, maxdepth = task$depth,
             minbucket = minb, subsample = 0.80,
             lambda = if (task$family == "gaussian") 1.0 else 2.0, nfolds = 3L)
}

ensemble_info <- function(fit, forest) {
  list(rho12 = fit$moments$cor[1, 2], w_alpha = fit$w_alpha, kl = fit$kl,
       iter = fit$iter, B = length(forest$trees), fallback = fit$fallback,
       q_binds = sum(fit$cuts$q_binds_lo | fit$cuts$q_binds_hi))
}

run_one <- function(task) {
  set.seed(task$seed)
  p_extra <- if (task$scenario == "p12") 7L else 0L
  X <- draw_x(task$n, task$rho, p_extra)
  y <- draw_y(X, task$family, task$scenario, task$b5)
  fam <- task$family
  M <- if (p_extra > 0) 4000L else 2200L
  minb <- if (task$study == "S1") max(8L, round(task$n / 60)) else 10L
  beta_full <- if (task$scenario == "linear" || task$scenario == "rare") {
    b <- if (task$scenario == "rare") c(-2.9, BETA[[fam]][2:6]) else c(BETA[[fam]][1:5], task$b5)
    setNames(b, c("(Intercept)", colnames(X)))
  } else true_projection(fam, task$scenario, task$rho, p_extra)
  out <- list()

  oracle <- ipd_glm(X, y, fam)
  gbt <- gbt_for(task, X, y, minb)
  fit <- tree_ipf(gbt, M = M, alpha = task$alpha, effects = task$study %in% c("S1", "S4", "S5"))
  info <- ensemble_info(fit, gbt)
  theta_I <- if (task$scenario %in% c("linear", "rare")) beta_full else
    true_projection(fam, task$scenario, task$rho, p_extra, fit$cuts)

  add <- function(method, est, se, target = theta_I, extra = list())
    out[[length(out) + 1L]] <<- pack(task, method, est, se, target, extra)

  add("Oracle", oracle$est, oracle$se, beta_full)
  rows_int <- in_box(X, fit$cuts)
  oi <- ipd_glm(X, y, fam, rows_int)
  add("Oracle-interior", oi$est, oi$se)
  add("Ensemble on true X", fhat_on_true_x(gbt, X, fit$cuts), rep(NA_real_, length(oracle$est)))
  add("TreeIPF", fit$est, fit$se, extra = info)
  # same TreeIPF estimate scored against the full-support target
  if (task$scenario == "nonlinear") add("TreeIPF vs full target", fit$est, fit$se, beta_full)

  if (task$study == "S3") {
    for (a in c(0, 0.02, 0.05, 0.10, 0.20)) {
      pr <- project_design(gbt, fit$Z, fit$q, fit$mu, fit$grid, alpha = a)
      add(sprintf("TreeIPF alpha=%.2f", a), pr$est, pr$se,
          extra = list(w_alpha = pr$w_alpha, q_binds = sum(pr$cuts$q_binds_lo | pr$cuts$q_binds_hi)))
    }
  }

  if (task$study %in% c("S1", "S4", "S5")) {
    ab <- ablation_grid(gbt, fit, alpha = task$alpha)
    for (k in names(ab)) add(paste0("A:", k), ab[[k]]$est, ab[[k]]$se)
    add("ALE", fit$effects$ale, fit$se)
    add("PD", fit$effects$pd, fit$se)
  }
  if (task$study == "S2") {
    d <- design_variant(gbt, fit, "marginal")
    pr <- project_design(gbt, d$Z, d$q, d$mu, fit$grid, task$alpha)
    add("Marg-IPF", pr$est, pr$se)
  }
  if (task$study == "S1") {
    bag <- fit_bagged_trees(X, y, fam, B = 100L, maxdepth = 5L, minbucket = minb, subsample = 0.80)
    fb <- tree_ipf(bag, M = M, alpha = task$alpha)
    add("TreeIPF-Bag", fb$est, fb$se, extra = ensemble_info(fb, bag))
  }
  # split counts and routed mass per covariate (identifiability diagnostics)
  ns <- fit$n_splits; om <- fit$routed_mass
  diag_dt <- data.table(study = task$study, family = fam, scenario = task$scenario,
                        n = task$n, rho = task$rho, depth = task$depth, b5 = task$b5,
                        rep = task$rep, term = gbt$var_names, n_splits = ns, routed_mass = om,
                        rec_mean = fit$moments$mean, rec_sd = fit$moments$sd,
                        emp_mean = colMeans(X), emp_sd = apply(X, 2, sd))
  list(est = rbindlist(out, fill = TRUE), diag = diag_dt)
}

# ------------------------------------------------------------------------------
# task lists (seeds are fixed per study, cell and replicate)
# ------------------------------------------------------------------------------

mk <- function(study, family, scenario = "linear", n, rho, depth = 3L, alpha = 0.05,
               b5 = 0, seed0) {
  g <- CJ(family = family, n = n, rho = rho, depth = depth, alpha = alpha, b5 = b5,
          scenario = scenario, rep = seq_len(R_REP))
  g[, study := study]
  g[, seed := seed0 + .I]
  g
}

tasks_for <- function(s) {
  fams <- c("gaussian", "binomial")
  switch(as.character(s),
    "1" = mk("S1", fams, n = c(300L, 600L, 1500L), rho = c(0, 0.4, 0.7), seed0 = 100000L),
    "2" = mk("S2", fams, n = 800L, rho = 0.7, depth = c(1L, 2L, 3L, 5L), seed0 = 500000L),
    "3" = mk("S3", fams, n = 800L, rho = 0.4, seed0 = 700000L),
    "4" = rbind(mk("S4", "gaussian", n = 600L, rho = 0.4, b5 = c(0, 0.05, 0.10, 0.15, 0.20, 0.30), seed0 = 800000L),
                mk("S4", "binomial", n = 600L, rho = 0.4, b5 = c(0, 0.075, 0.15, 0.225, 0.30, 0.45), seed0 = 810000L)),
    "5" = rbind(mk("S5", fams, "nonlinear", n = c(600L, 1500L), rho = 0.4, seed0 = 900000L),
                mk("S5", "binomial", "rare", n = 1500L, rho = 0.4, seed0 = 910000L),
                mk("S5", fams, "p12", n = 1500L, rho = 0.4, seed0 = 920000L)))
}

rep_file <- function(task) {
  file.path(root_dir, "results", "reps", task$study,
            sprintf("%s_%s_n%d_r%.1f_d%d_a%.2f_b%.3f_%03d.rds", task$family, task$scenario,
                    task$n, task$rho, task$depth, task$alpha, task$b5, task$rep))
}

run_task_safe <- function(i, tasks) {
  task <- as.list(tasks[i])
  f <- rep_file(task)
  if (file.exists(f)) return(f)
  res <- tryCatch(run_one(task), error = function(e) list(error = conditionMessage(e)))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  saveRDS(res, f)
  f
}

if (sys.nframe() == 0L) {
  for (s in studies) {
    tasks <- tasks_for(s)
    todo <- which(!file.exists(vapply(seq_len(nrow(tasks)), function(i) rep_file(as.list(tasks[i])), "")))
    cat(sprintf("[%s] study %d: %d tasks, %d to run, %d workers\n",
                format(Sys.time(), "%H:%M"), s, nrow(tasks), length(todo), workers))
    if (!length(todo)) next
    cl <- makeCluster(workers)
    clusterExport(cl, c("lib_file", "root_dir"))
    clusterEvalQ(cl, { suppressPackageStartupMessages(library(data.table)); source(lib_file); NULL })
    clusterExport(cl, c("BETA", "draw_x", "true_eta", "draw_y", "BIG", "big_sample",
                        "true_projection", "pack", "gbt_for", "ensemble_info", "run_one",
                        "rep_file", "run_task_safe", "tasks"))
    invisible(parLapplyLB(cl, todo, run_task_safe, tasks = tasks))
    stopCluster(cl)
    cat(sprintf("[%s] study %d done\n", format(Sys.time(), "%H:%M"), s))
  }
}
