# ==============================================================================
# 08_real_libraries.R  --  TreeIPF on model files saved by real libraries (Study 7)
#
# Same cohorts as Study 1, cell n = 1500, rho = 0.4 (same seeds), R = 100.
#   XGBoost  : JSON file (xgb.save), Hessian covers, no root deviance;
#              sigma^2 from the cross-validated RMSE stored as a model attribute
#   sklearn  : RandomForest, max_features = 0.5 (mtry < p), bootstrap, JSON export
#   ranger   : node.stats = TRUE (node counts kept), mtry = 2, sigma^2 from OOB error
#
# usage: Rscript 08_real_libraries.R [workers]
# ==============================================================================

suppressPackageStartupMessages({ library(parallel); library(data.table) })
root_dir <- normalizePath(Sys.getenv("TREEIPF_ROOT", "."), winslash = "/")
sim_file <- file.path(root_dir, "code", "01_run_simulations.R")
source(sim_file)
args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args) >= 1) as.integer(args[1]) else 8L
R_LIB <- 100L

lib_dir <- file.path(root_dir, "results", "libmodels")
csv_dir <- file.path(lib_dir, "csv"); json_dir <- file.path(lib_dir, "sklearn")
xgb_dir <- file.path(lib_dir, "xgboost"); rgr_dir <- file.path(lib_dir, "ranger")
for (d in c(csv_dir, json_dir, xgb_dir, rgr_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

tasks_l <- tasks_for(1)[n == 1500L & rho == 0.4 & rep <= R_LIB]
tag <- function(task) sprintf("%s_%03d", task$family, task$rep)

cohort <- function(task) {
  set.seed(task$seed)
  X <- draw_x(task$n, task$rho, 0L)
  y <- draw_y(X, task$family, "linear", 0)
  list(X = X, y = y)
}

fit_xgb <- function(X, y, family, path) {
  library(xgboost)
  obj <- if (family == "gaussian") "reg:squarederror" else "binary:logistic"
  par <- list(objective = obj, eta = 0.1, max_depth = 3, subsample = 0.8,
              lambda = 1, min_child_weight = if (family == "gaussian") 10 else 2, nthread = 1)
  dm <- xgb.DMatrix(X, label = y)
  cv <- xgb.cv(par, dm, nrounds = 400, nfold = 5, early_stopping_rounds = 20, verbose = 0)
  best <- cv$early_stop$best_iteration
  if (is.null(best)) best <- cv$best_iteration
  ev <- as.data.frame(cv$evaluation_log)
  bst <- xgb.train(par, dm, nrounds = best, verbose = 0)
  if (family == "gaussian") xgb.attr(bst, "best_score") <- ev$test_rmse_mean[best]
  xgb.save(bst, path)
  best
}

fit_ranger <- function(X, y, family, path) {
  d <- data.frame(y = if (family == "binomial") factor(y) else y, X)
  rf <- ranger::ranger(y ~ ., d, num.trees = 200, mtry = 2, min.node.size = 10,
                       probability = family == "binomial", node.stats = TRUE,
                       num.threads = 1, seed = 1)
  saveRDS(rf, path)
}

lib_one <- function(i, tasks) {
  task <- as.list(tasks[i])
  f <- file.path(root_dir, "results", "reps", "S7", paste0(tag(task), ".rds"))
  if (file.exists(f)) return(f)
  res <- tryCatch({
    co <- cohort(task); X <- co$X; y <- co$y; fam <- task$family
    target <- setNames(c(BETA[[fam]][1:5], 0), c("(Intercept)", colnames(X)))
    out <- list()
    add <- function(method, est, se, extra = list())
      out[[length(out) + 1L]] <<- pack(task, method, est, se, target, extra)
    o <- ipd_glm(X, y, fam); add("Oracle", o$est, o$se)

    xp <- file.path(xgb_dir, paste0(tag(task), ".json"))
    if (!file.exists(xp)) fit_xgb(X, y, fam, xp)
    rp <- file.path(rgr_dir, paste0(tag(task), ".rds"))
    if (!file.exists(rp)) fit_ranger(X, y, fam, rp)
    sp <- file.path(json_dir, paste0(tag(task), ".json"))

    forests <- list(
      XGBoost = read_xgb_json(xp, subsample = 0.8),
      sklearn = read_sklearn_json(sp),
      ranger  = {
        rf <- readRDS(rp); fr <- read_ranger(rf, colnames(X))
        if (fam == "gaussian") fr$sigma2 <- rf$prediction.error
        fr
      })
    for (nm in names(forests)) {
      fo <- forests[[nm]]
      Mg <- if (fo$mode == "bagged") 5000L else 2200L
      modes <- if (nm == "XGBoost" && fam == "binomial") c("hessian", "count") else "hessian"
      for (cm in modes) {
        set.seed(task$seed + 77L)
        ft <- tree_ipf(fo, M = Mg, cover_mode = cm)
        lab <- if (cm == "count") paste(nm, "(covers as counts)") else nm
        add(paste("TreeIPF", lab), ft$est, ft$se,
            extra = list(n_hat = fo$n_hat, sigma2 = ft$sigma2, rho12 = ft$moments$cor[1, 2],
                         miss_mass = ft$miss_mass, B = length(fo$trees)))
        if (cm == "hessian") {
          d <- design_variant(fo, ft, "uniform")
          pr <- project_design(fo, d$Z, d$q, d$mu, ft$grid, 0.05)
          add(paste("Unif-rootn", nm), pr$est, pr$se)
        }
      }
    }
    rbindlist(out, fill = TRUE)
  }, error = function(e) list(error = conditionMessage(e)))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  saveRDS(res, f)
  f
}

if (sys.nframe() == 0L) {
  # 1) cohorts to csv for the python exporter
  for (i in seq_len(nrow(tasks_l))) {
    task <- as.list(tasks_l[i]); cf <- file.path(csv_dir, paste0(tag(task), ".csv"))
    if (!file.exists(cf)) { co <- cohort(task); fwrite(data.frame(co$X, y = co$y), cf) }
  }
  # 2) scikit-learn forests
  system2("python", c(shQuote(file.path(root_dir, "code", "sklearn_export.py")), shQuote(csv_dir), shQuote(json_dir)))
  # 3) xgboost / ranger fits and TreeIPF on all three files
  todo <- seq_len(nrow(tasks_l))
  cl <- makeCluster(workers)
  clusterExport(cl, c("sim_file", "root_dir"))
  clusterEvalQ(cl, { suppressPackageStartupMessages(library(data.table)); source(sim_file); NULL })
  clusterExport(cl, c("lib_dir", "csv_dir", "json_dir", "xgb_dir", "rgr_dir", "tag", "cohort",
                      "fit_xgb", "fit_ranger", "lib_one", "tasks_l"))
  invisible(parLapplyLB(cl, todo, lib_one, tasks = tasks_l))
  stopCluster(cl)
  cat("real-library study done\n")
}
