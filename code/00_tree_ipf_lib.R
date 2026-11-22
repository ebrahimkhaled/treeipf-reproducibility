# ==============================================================================
# 00_tree_ipf_lib.R  --  TreeIPF core library (version 3)
#
# Everything works on one plain tree format, so the same code reads
#   rpart ensembles (our own booster / bagger), XGBoost JSON files,
#   scikit-learn forests (exported to JSON by sklearn_export.py) and ranger.
#
# A tree is a list of node arrays (1-based):
#   feat  : feature index (0 = leaf)       thr   : split threshold
#   left  : row of the left child          right : row of the right child
#   value : leaf value (prediction or boosting step)
#   size  : training count (or Hessian cover) in the node
#   le    : TRUE if "x <= thr" goes left, FALSE if "x < thr" goes left
# ==============================================================================

local({
  libs <- strsplit(Sys.getenv("TREEIPF_RLIB", ""), ";")[[1]]   # optional extra R library folders
  libs <- libs[dir.exists(libs)]
  if (length(libs)) .libPaths(c(libs, .libPaths()))
})
suppressPackageStartupMessages({
  library(rpart)
  library(MASS)
})

sig   <- function(u) 1 / (1 + exp(-pmin(pmax(u, -30), 30)))
clamp <- function(p, e = 1e-5) pmin(pmax(p, e), 1 - e)

# ------------------------------------------------------------------------------
# 1. tree format, routing, prediction
# ------------------------------------------------------------------------------

new_forest <- function(trees, var_names, family, mode, base_score, lr = 1,
                       n_hat, tss = NA_real_, sigma2 = NA_real_,
                       size_type = "count", source = "rpart", hyper = list()) {
  structure(list(trees = trees, var_names = var_names, family = family,
                 mode = mode, base_score = base_score, lr = lr,
                 n_hat = n_hat, tss = tss, sigma2 = sigma2,
                 size_type = size_type, source = source, hyper = hyper),
            class = "tree_forest")
}

# leaf row reached by every row of Z
route_tree <- function(tr, Z) {
  node <- rep(1L, nrow(Z))
  todo <- which(tr$feat[node] > 0L)
  while (length(todo)) {
    nd <- node[todo]
    x  <- Z[cbind(todo, tr$feat[nd])]
    go_left <- if (tr$le) x <= tr$thr[nd] else x < tr$thr[nd]
    node[todo] <- ifelse(go_left, tr$left[nd], tr$right[nd])
    todo <- todo[tr$feat[node[todo]] > 0L]
  }
  node
}

# raw per-tree leaf values on Z (M x B) plus leaf rows
forest_leaves <- function(forest, Z) {
  B <- length(forest$trees)
  leaf <- matrix(0L, nrow(Z), B)
  vals <- matrix(0, nrow(Z), B)
  for (b in seq_len(B)) {
    tr <- forest$trees[[b]]
    leaf[, b] <- route_tree(tr, Z)
    vals[, b] <- tr$value[leaf[, b]]
  }
  list(leaf = leaf, vals = vals)
}

# prediction on the link scale (identity for gaussian, logit for binomial)
forest_link <- function(forest, Z, vals = NULL) {
  if (is.null(vals)) vals <- forest_leaves(forest, Z)$vals
  if (forest$mode == "boosted") {
    forest$base_score + forest$lr * rowSums(vals)
  } else {
    m <- rowMeans(vals)
    if (forest$family == "binomial") qlogis(clamp(m)) else m
  }
}

forest_mean <- function(forest, Z, vals = NULL) {
  eta <- forest_link(forest, Z, vals)
  if (forest$family == "binomial") sig(eta) else eta
}

# ------------------------------------------------------------------------------
# 2. readers for the four sources
# ------------------------------------------------------------------------------

rpart_to_tree <- function(tr, var_names) {
  fr <- tr$frame
  ids <- as.integer(rownames(fr))
  leaf <- fr$var == "<leaf>"
  K <- nrow(fr)
  feat  <- ifelse(leaf, 0L, match(as.character(fr$var), var_names))
  thr   <- rep(NA_real_, K)
  left  <- rep(NA_integer_, K)
  right <- rep(NA_integer_, K)
  inner <- which(!leaf)
  sp <- tr$splits
  for (k in seq_along(inner)) {
    r <- inner[k]
    thr[r] <- sp[k, "index"]
    lo <- match(2L * ids[r], ids); hi <- match(2L * ids[r] + 1L, ids)
    # rpart: ncat = -1 means "x < index" goes to the left child (2k)
    if (sp[k, "ncat"] < 0) { left[r] <- lo; right[r] <- hi } else { left[r] <- hi; right[r] <- lo }
  }
  list(feat = as.integer(feat), thr = thr, left = left, right = right,
       value = fr$yval, size = fr$n, le = FALSE)
}

# XGBoost model saved with xgb.save(model, "file.json") / save_model("file.json")
read_xgb_json <- function(path, subsample = 1, sigma2 = NA_real_, var_names = NULL) {
  j <- jsonlite::fromJSON(path, simplifyVector = TRUE, simplifyDataFrame = FALSE)
  L <- j$learner
  obj <- L$objective$name
  family <- if (grepl("logistic", obj)) "binomial" else "gaussian"
  bs <- as.numeric(gsub("\\[|\\]", "", L$learner_model_param$base_score))
  base_margin <- if (family == "binomial") qlogis(bs) else bs
  fn <- unlist(L$feature_names)
  if (is.null(var_names)) var_names <- fn
  trees <- lapply(L$gradient_booster$model$trees, function(t) {
    lc <- unlist(t$left_children); rc <- unlist(t$right_children)
    is_leaf <- lc < 0
    feat <- ifelse(is_leaf, 0L, match(fn[unlist(t$split_indices) + 1L], var_names))
    list(feat = as.integer(feat), thr = unlist(t$split_conditions),
         left = ifelse(is_leaf, NA_integer_, lc + 1L),
         right = ifelse(is_leaf, NA_integer_, rc + 1L),
         value = ifelse(is_leaf, unlist(t$split_conditions), 0),
         size = unlist(t$sum_hessian), le = FALSE)
  })
  # cohort size from the first tree's root cover (Remark 1):
  #   squared error: hessian = 1 per row;  logistic: hessian = p0 (1 - p0) per row
  h0 <- if (family == "binomial") bs * (1 - bs) else 1
  n_hat <- round(trees[[1]]$size[1] / (subsample * h0))
  attrs <- L$attributes
  if (is.na(sigma2) && family == "gaussian" && !is.null(attrs$best_score))
    sigma2 <- as.numeric(attrs$best_score)^2
  new_forest(trees, var_names, family, "boosted", base_margin, lr = 1,
             n_hat = n_hat, tss = NA_real_, sigma2 = sigma2,
             size_type = if (family == "binomial") "cover" else "count",
             source = "xgboost", hyper = list(subsample = subsample))
}

# scikit-learn forest exported by code/sklearn_export.py
read_sklearn_json <- function(path) {
  j <- jsonlite::fromJSON(path, simplifyVector = TRUE, simplifyDataFrame = FALSE)
  family <- if (j$task == "classification") "binomial" else "gaussian"
  vn <- unlist(j$feature_names)
  trees <- lapply(j$trees, function(t) {
    lc <- unlist(t$left); rc <- unlist(t$right)
    is_leaf <- lc < 0
    list(feat = as.integer(ifelse(is_leaf, 0L, unlist(t$feature) + 1L)),
         thr = unlist(t$threshold),
         left = ifelse(is_leaf, NA_integer_, lc + 1L),
         right = ifelse(is_leaf, NA_integer_, rc + 1L),
         value = unlist(t$value), size = unlist(t$wn), le = TRUE,
         root_impurity = t$impurity[[1]])
  })
  # bootstrap draws n rows per tree, so the weighted root count is n itself
  frac <- if (is.null(j$max_samples)) 1 else j$max_samples
  n_hat <- round(trees[[1]]$size[1] / frac)
  tss <- if (family == "gaussian") {
    mean(sapply(trees, function(t) t$root_impurity * t$size[1])) / frac
  } else NA_real_
  new_forest(trees, vn, family, "bagged", NA_real_, 1, n_hat = n_hat, tss = tss,
             source = "sklearn", hyper = list(max_features = j$max_features))
}

# ranger forest grown with node.stats = TRUE (otherwise node counts are not kept)
read_ranger <- function(rf, var_names) {
  if (is.null(rf$forest$num.samples.nodes))
    stop("ranger forest has no node counts; refit with node.stats = TRUE")
  family <- if (rf$treetype == "Regression") "gaussian" else "binomial"
  trees <- lapply(seq_len(rf$num.trees), function(b) {
    ti <- ranger::treeInfo(rf, b)
    pcol <- if (family == "binomial") grep("^pred\\.", names(ti), value = TRUE) else "prediction"
    val <- if (family == "binomial") ti[[pcol[length(pcol)]]] else ti$prediction
    leaf <- ti$terminal
    list(feat = as.integer(ifelse(leaf, 0L, match(ti$splitvarName, var_names))),
         thr = ti$splitval,
         left = ifelse(leaf, NA_integer_, ti$leftChild + 1L),
         right = ifelse(leaf, NA_integer_, ti$rightChild + 1L),
         value = val, size = ti$numSamples, le = TRUE)
  })
  frac <- rf$call$sample.fraction
  frac <- if (is.null(frac)) 1 else eval(frac)
  new_forest(trees, var_names, family, "bagged", NA_real_, 1,
             n_hat = round(trees[[1]]$size[1] / frac),
             source = "ranger", hyper = list(mtry = rf$mtry))
}

# ------------------------------------------------------------------------------
# 3. our own ensembles (rpart trees, data thrown away on return)
# ------------------------------------------------------------------------------

grow_rpart <- function(g, X_df, minbucket, maxdepth) {
  rpart(g ~ ., data = data.frame(g = g, X_df),
        control = rpart.control(cp = 0, minbucket = minbucket, maxdepth = maxdepth,
                                maxcompete = 0, maxsurrogate = 0, xval = 0))
}

# one Newton-boosting run with a fixed number of rounds; returns rpart trees
boost_rounds <- function(X_df, y, family, B, lr, maxdepth, minbucket, subsample,
                         lambda, X_pred = NULL) {
  n <- nrow(X_df)
  n_sub <- max(20L, floor(subsample * n))
  rm_all <- rpart:::rpart.matrix(X_df)
  if (family == "gaussian") base <- mean(y) else base <- qlogis(min(max(mean(y), 1e-4), 1 - 1e-4))
  fc <- rep(base, n)
  pred_path <- if (!is.null(X_pred)) matrix(0, nrow(X_pred), B) else NULL
  fp <- if (!is.null(X_pred)) rep(base, nrow(X_pred)) else NULL
  rm_pred <- if (!is.null(X_pred)) rpart:::rpart.matrix(X_pred) else NULL
  trees <- vector("list", B)
  for (b in seq_len(B)) {
    si <- if (subsample < 1) sample.int(n, n_sub) else seq_len(n)
    if (family == "gaussian") {
      g <- y[si] - fc[si]; h <- rep(1, length(si))
    } else {
      pc <- sig(fc[si]); g <- y[si] - pc; h <- pc * (1 - pc)
    }
    tr <- grow_rpart(g, X_df[si, , drop = FALSE], minbucket, maxdepth)
    w <- rpart:::pred.rpart(tr, rpart:::rpart.matrix(X_df[si, , drop = FALSE]))
    sg <- tapply(g, w, sum); sh <- tapply(h, w, sum)
    nd <- as.integer(names(sg))
    tr$frame$yval[nd] <- sg / (sh + lambda)
    fc <- fc + lr * tr$frame$yval[rpart:::pred.rpart(tr, rm_all)]
    if (!is.null(X_pred)) {
      fp <- fp + lr * tr$frame$yval[rpart:::pred.rpart(tr, rm_pred)]
      pred_path[, b] <- fp
    }
    trees[[b]] <- tr
  }
  list(trees = trees, base = base, pred_path = pred_path, n_sub = n_sub)
}

fit_cv_gbt <- function(X, y, family = c("gaussian", "binomial"),
                       B_max = 150L, lr = 0.06, maxdepth = 3L,
                       minbucket = 10L, subsample = 0.80, lambda = 1.5,
                       nfolds = 5L, fixed_B = NULL, min_B = 15L) {
  family <- match.arg(family)
  X_df <- as.data.frame(X)
  n <- nrow(X_df)
  if (is.null(fixed_B)) {
    folds <- sample(rep(seq_len(nfolds), length.out = n))
    oof <- matrix(0, n, B_max)
    for (k in seq_len(nfolds)) {
      tr_i <- which(folds != k); va_i <- which(folds == k)
      run <- boost_rounds(X_df[tr_i, , drop = FALSE], y[tr_i], family, B_max, lr,
                          maxdepth, minbucket, 1, lambda, X_pred = X_df[va_i, , drop = FALSE])
      oof[va_i, ] <- run$pred_path
    }
    loss <- if (family == "gaussian") colMeans((y - oof)^2) else {
      pv <- sig(oof); colMeans(-y * log(pv + 1e-12) - (1 - y) * log(1 - pv + 1e-12))
    }
    B_opt <- max(min_B, which.min(loss))
  } else B_opt <- as.integer(fixed_B)

  run <- boost_rounds(X_df, y, family, B_opt, lr, maxdepth, minbucket, subsample, lambda)
  vn <- colnames(X_df)
  trees <- lapply(run$trees, rpart_to_tree, var_names = vn)
  # root of tree 1 is fitted to y - mean(y): its deviance is the subsample TSS
  tss <- if (family == "gaussian") run$trees[[1]]$frame$dev[1] / subsample else NA_real_
  new_forest(trees, vn, family, "boosted", run$base, lr,
             n_hat = round(trees[[1]]$size[1] / subsample), tss = tss,
             hyper = list(type = "gbt", B = B_opt, lr = lr, maxdepth = maxdepth,
                          minbucket = minbucket, subsample = subsample, lambda = lambda))
}

fit_bagged_trees <- function(X, y, family = c("gaussian", "binomial"),
                             B = 100L, maxdepth = 5L, minbucket = 8L, subsample = 0.80) {
  family <- match.arg(family)
  X_df <- as.data.frame(X)
  n <- nrow(X_df)
  n_sub <- max(20L, floor(subsample * n))
  vn <- colnames(X_df)
  raw <- lapply(seq_len(B), function(b) {
    si <- sample.int(n, n_sub)
    grow_rpart(y[si], X_df[si, , drop = FALSE], minbucket, maxdepth)
  })
  trees <- lapply(raw, rpart_to_tree, var_names = vn)
  tss <- if (family == "gaussian") mean(sapply(raw, function(t) t$frame$dev[1])) / subsample else NA_real_
  new_forest(trees, vn, family, "bagged", NA_real_, 1,
             n_hat = round(trees[[1]]$size[1] / subsample), tss = tss,
             hyper = list(type = "bag", B = B, maxdepth = maxdepth,
                          minbucket = minbucket, subsample = subsample))
}

refit_like <- function(forest, X, y) {
  h <- forest$hyper
  if (identical(h$type, "gbt")) {
    fit_cv_gbt(X, y, forest$family, lr = h$lr, maxdepth = h$maxdepth, minbucket = h$minbucket,
               subsample = h$subsample, lambda = h$lambda, fixed_B = h$B)
  } else if (identical(h$type, "bag")) {
    fit_bagged_trees(X, y, forest$family, B = h$B, maxdepth = h$maxdepth,
                     minbucket = h$minbucket, subsample = h$subsample)
  } else stop("no refit rule for this forest")
}

# ------------------------------------------------------------------------------
# 4. what the trees reveal: split sets, variable types, identifiability
# ------------------------------------------------------------------------------

split_sets <- function(forest) {
  p <- length(forest$var_names)
  sp <- vector("list", p); cnt <- integer(p)
  for (tr in forest$trees) {
    inner <- which(tr$feat > 0L)
    for (r in inner) {
      j <- tr$feat[r]
      sp[[j]] <- c(sp[[j]], tr$thr[r]); cnt[j] <- cnt[j] + 1L
    }
  }
  list(sets = lapply(sp, function(s) sort(unique(s))), n_splits = cnt)
}

# rules (documented in the paper, Methods / Additional file):
#   no split            -> unsplit (dropped from the design, estimate NA)
#   one split at 0.5     -> binary {0, 1}
#   all splits at k+0.5  -> integer score on floor(min)..ceiling(max)  (span <= 250)
#   otherwise            -> continuous
# XGBoost splits AT observed values ("x < v"), not halfway between them, so there
# a binary 0/1 predictor shows one split at 1 and an integer score splits at integers
at_values <- function(forest) identical(forest$source, "xgboost")

var_types <- function(sets, at_val = FALSE) {
  vapply(sets, function(s) {
    if (!length(s)) return("unsplit")
    if (length(s) == 1L && (abs(s - 0.5) < 1e-5 || (at_val && abs(s - 1) < 1e-5))) return("binary")
    off <- if (at_val) 0 else 0.5
    if (all(abs((s - off) - round(s - off)) < 1e-5) && diff(range(s)) <= 250) return("integer")
    "continuous"
  }, character(1))
}

# share of training mass whose path uses variable j, averaged over trees
split_routed_mass <- function(forest) {
  p <- length(forest$var_names)
  om <- numeric(p)
  for (tr in forest$trees) {
    K <- length(tr$feat)
    on <- matrix(FALSE, K, p)
    for (r in seq_len(K)) if (tr$feat[r] > 0L) {
      fl <- on[r, ]; fl[tr$feat[r]] <- TRUE
      on[tr$left[r], ] <- fl; on[tr$right[r], ] <- fl
    }
    lv <- which(tr$feat == 0L)
    om <- om + colSums(on[lv, , drop = FALSE] * (tr$size[lv] / tr$size[1]))
  }
  setNames(om / length(forest$trees), forest$var_names)
}

# ------------------------------------------------------------------------------
# 5. step 1: grid and multiplicative I-projection (SMART)
# ------------------------------------------------------------------------------

build_grid <- function(forest, M = 2200L, sets = NULL, design = c("strat", "uniform")) {
  design <- match.arg(design)
  if (is.null(sets)) sets <- split_sets(forest)$sets
  vt <- var_types(sets, at_values(forest))
  p <- length(sets)
  Z <- matrix(0, M, p, dimnames = list(NULL, forest$var_names))
  lo <- hi <- numeric(p)
  for (j in seq_len(p)) {
    s <- sets[[j]]
    if (vt[j] == "unsplit") {
      Z[, j] <- 0; lo[j] <- hi[j] <- 0           # never used by any tree
    } else if (vt[j] == "binary") {
      Z[, j] <- sample(c(0, 1), M, TRUE); lo[j] <- 0; hi[j] <- 1
    } else if (vt[j] == "integer") {
      v <- if (at_values(forest)) seq(round(min(s)) - 1, round(max(s))) else seq(floor(min(s)), ceiling(max(s)))
      Z[, j] <- v[sample.int(length(v), M, TRUE)]; lo[j] <- min(v); hi[j] <- max(v)
    } else {
      if (length(s) >= 4L) {
        iqr <- diff(quantile(s, c(0.25, 0.75), names = FALSE))
        if (iqr <= 0) iqr <- max(diff(range(s)) * 0.5, 0.5)
        pad <- getOption("treeipf.pad", 0.25) * iqr      # tail rule; varied in 12_checks.R
      } else pad <- 2 * getOption("treeipf.pad", 0.25) * max(diff(range(s)), 1)
      knots <- c(s[1] - pad, s, s[length(s)] + pad)
      lo[j] <- knots[1]; hi[j] <- knots[length(knots)]
      if (design == "strat") {
        k <- sample.int(length(knots) - 1L, M, TRUE)
        Z[, j] <- knots[k] + runif(M) * (knots[k + 1L] - knots[k])
      } else Z[, j] <- runif(M, lo[j], hi[j])
    }
  }
  list(Z = Z, sets = sets, var_type = vt, lo = lo, hi = hi)
}

# Hessian weights per grid point and tree for boosted logistic covers (Remark 1)
cover_weights <- function(forest, vals) {
  if (forest$size_type != "cover") return(NULL)
  eta <- forest$base_score + forest$lr * cbind(0, t(apply(vals, 1, cumsum)))[, seq_len(ncol(vals)), drop = FALSE]
  p <- sig(eta); p * (1 - p)
}

smart_fit <- function(forest, leaves, H = NULL, gamma = 1, max_iter = 150L, tol = 1e-6) {
  M <- nrow(leaves$leaf); B <- ncol(leaves$leaf)
  offs <- c(0L, cumsum(vapply(forest$trees, function(t) length(t$feat), 1L)))
  gid <- as.vector(sweep(leaves$leaf, 2, offs[-length(offs)], "+"))
  target <- unlist(lapply(forest$trees, function(t) t$size / t$size[1]))
  nl <- length(target)
  tg <- target[gid]
  hit <- unique(gid)
  lv_all <- unlist(lapply(seq_len(B), function(b) offs[b] + which(forest$trees[[b]]$feat == 0L)))
  miss <- setdiff(lv_all, hit)
  # sparse leaf-by-point incidence (Hessian-weighted for logistic covers)
  hv <- if (is.null(H)) rep(1, M * B) else as.vector(H)
  A <- Matrix::sparseMatrix(i = gid, j = rep(seq_len(M), B), x = hv, dims = c(nl, M))
  tree_of <- findInterval(seq_len(nl) - 1, offs[-1]) + 1L
  q <- rep(1 / M, M)
  kl_old <- Inf
  for (it in seq_len(max_iter)) {
    mass <- as.vector(A %*% q)
    if (!is.null(H)) mass <- mass / rowsum(mass, tree_of, reorder = TRUE)[tree_of, 1]
    kl <- sum(mass[hit] * log((mass[hit] + 1e-12) / (target[hit] + 1e-12))) / B
    lr <- log((target + 1e-9) / (mass + 1e-9))
    lm <- as.vector(Matrix::crossprod(A != 0, lr)) / B
    q <- q * exp(gamma * lm); q <- q / sum(q)
    if (abs(kl_old - kl) < tol * max(kl, 1e-4)) break
    kl_old <- kl
  }
  list(q = q, kl = kl, iter = it, miss_leaves = length(miss), miss_mass = sum(target[miss]) / B)
}

# ------------------------------------------------------------------------------
# 6. steps 2-3: interior projection and variance
# ------------------------------------------------------------------------------

interior_cuts <- function(Z, q, sets, vt, alpha) {
  p <- ncol(Z)
  cuts <- data.frame(var = colnames(Z), lo = -Inf, hi = Inf, trimmed = FALSE,
                     q_binds_lo = FALSE, q_binds_hi = FALSE)
  if (alpha <= 0) return(cuts)
  a_j <- alpha / max(1L, sum(vt == "continuous"))
  for (j in seq_len(p)) {
    s <- sets[[j]]
    if (vt[j] == "continuous" && length(s) >= 6L) {
      o <- order(Z[, j]); cw <- cumsum(q[o])
      ql <- Z[o[which.max(cw >= a_j)], j]; qh <- Z[o[which.max(cw >= 1 - a_j)], j]
      lo <- max(s[1], ql); hi <- min(s[length(s)], qh)
      if (lo < hi) {
        cuts$lo[j] <- lo; cuts$hi[j] <- hi; cuts$trimmed[j] <- TRUE
        cuts$q_binds_lo[j] <- ql > s[1]; cuts$q_binds_hi[j] <- qh < s[length(s)]
      }
    }
  }
  cuts
}

in_box <- function(X, cuts) {
  ok <- rep(TRUE, nrow(X))
  for (j in which(cuts$trimmed)) ok <- ok & X[, j] >= cuts$lo[j] & X[, j] <= cuts$hi[j]
  ok
}

wls_glm <- function(D, w, mu, family, iters = 30L) {
  if (family == "gaussian") return(solve(crossprod(D, w * D), crossprod(D, w * mu))[, 1])
  b <- c(qlogis(min(max(sum(w * mu) / sum(w), 1e-4), 1 - 1e-4)), rep(0, ncol(D) - 1L))
  for (it in seq_len(iters)) {
    pm <- sig(as.vector(D %*% b))
    step <- solve(crossprod(D, (w * pm * (1 - pm)) * D) + diag(1e-10, ncol(D)),
                  crossprod(D, w * (mu - pm)))[, 1]
    b <- b + step
    if (max(abs(step)) < 1e-9) break
  }
  b
}

# project a fitted mean surface onto the working GLM over a weighted design
#   se_mode "rootn": n from the root node, sigma^2 from root TSS (TreeIPF step 3)
#   se_mode "grid" : treats the M design points as the sample (naive distillation)
project_design <- function(forest, Z, q, mu, info, alpha = 0.05, trim = TRUE,
                           se_mode = c("rootn", "grid"), q_full = NULL, Z_full = NULL) {
  se_mode <- match.arg(se_mode)
  vt <- info$var_type
  keep <- which(vt != "unsplit")
  cuts <- interior_cuts(Z, q, info$sets, vt, if (trim) alpha else 0)
  inside <- in_box(Z, cuts)
  fallback <- sum(inside) < 5L * (length(keep) + 1L)
  if (fallback) inside <- rep(TRUE, nrow(Z))
  w_a <- sum(q[inside])
  qi <- q[inside] / w_a
  D <- cbind(1, Z[inside, keep, drop = FALSE])
  b <- wls_glm(D, qi, mu[inside], forest$family)

  n_hat <- forest$n_hat
  if (se_mode == "rootn") {
    n_int <- n_hat * w_a
    if (forest$family == "gaussian") {
      mz <- colSums(q * Z[, keep, drop = FALSE])
      Zc <- sweep(Z[, keep, drop = FALSE], 2, mz)
      S <- crossprod(Zc, q * Zc)
      if (!is.na(forest$tss)) {
        s2 <- (forest$tss - n_hat * drop(t(b[-1]) %*% S %*% b[-1])) / max(n_hat - length(keep) - 1, 10)
        # recovered signal larger than the total: fall back to the total variance (an upper bound)
        if (s2 <= 0) s2 <- forest$tss / (n_hat - 1)
      } else s2 <- forest$sigma2
      V <- s2 * solve(n_int * crossprod(D, qi * D))
    } else {
      s2 <- NA_real_
      pm <- sig(as.vector(D %*% b))
      V <- solve(n_int * crossprod(D, (qi * pm * (1 - pm)) * D))
    }
  } else {
    Mi <- sum(inside); wi <- qi * Mi
    if (forest$family == "gaussian") {
      r <- mu[inside] - as.vector(D %*% b)
      s2 <- sum(wi * r^2) / max(Mi - ncol(D), 1)
      V <- s2 * solve(crossprod(D, wi * D))
    } else {
      s2 <- NA_real_
      pm <- sig(as.vector(D %*% b))
      V <- solve(crossprod(D, (wi * pm * (1 - pm)) * D))
    }
  }
  p <- ncol(Z)
  est <- se <- rep(NA_real_, p + 1L)
  est[c(1L, keep + 1L)] <- b
  se[keep + 1L] <- sqrt(pmax(diag(V)[-1], 0))
  # intercept re-anchored over the full design to the root mean; reported without SE
  Zf <- if (is.null(Z_full)) Z else Z_full
  qf <- if (is.null(q_full)) q else q_full
  root_mean <- if (forest$family == "gaussian") {
    if (forest$mode == "boosted") forest$base_score else sum(qf * forest_mean(forest, Zf))
  } else {
    if (forest$mode == "boosted" && forest$source != "xgboost") sig(forest$base_score) else sum(qf * forest_mean(forest, Zf))
  }
  lp <- as.vector(Zf[, keep, drop = FALSE] %*% b[-1])
  est[1] <- if (forest$family == "gaussian") root_mean - sum(qf * lp) else
    tryCatch(uniroot(function(a) sum(qf * sig(a + lp)) - root_mean, c(-50, 50))$root,
             error = function(e) b[1])
  names(est) <- names(se) <- c("(Intercept)", forest$var_names)
  list(est = est, se = se, cuts = cuts, w_alpha = w_a, fallback = fallback, sigma2 = s2)
}

weighted_moments <- function(Z, q) {
  mz <- colSums(q * Z); Zc <- sweep(Z, 2, mz)
  S <- crossprod(Zc, q * Zc); sd <- sqrt(pmax(diag(S), 1e-12))
  list(mean = mz, cov = S, sd = sd, cor = S / outer(sd, sd))
}

# ------------------------------------------------------------------------------
# 7. the full TreeIPF fit (and every variant used in the ablation)
# ------------------------------------------------------------------------------

tree_ipf <- function(forest, M = 2200L, alpha = 0.05, gamma = 1, max_iter = 150L,
                     cover_mode = c("hessian", "count"), effects = FALSE) {
  cover_mode <- match.arg(cover_mode)
  ss <- split_sets(forest)
  gr <- build_grid(forest, M, ss$sets)
  lv <- forest_leaves(forest, gr$Z)
  H <- if (cover_mode == "hessian") cover_weights(forest, lv$vals) else NULL
  sm <- smart_fit(forest, lv, H, gamma, max_iter)
  mu <- forest_mean(forest, gr$Z, lv$vals)
  fit <- project_design(forest, gr$Z, sm$q, mu, gr, alpha)
  mo <- weighted_moments(gr$Z, sm$q)
  out <- c(fit, list(Z = gr$Z, q = sm$q, mu = mu, grid = gr, moments = mo,
                     kl = sm$kl, iter = sm$iter, miss_leaves = sm$miss_leaves,
                     miss_mass = sm$miss_mass, n_splits = ss$n_splits,
                     routed_mass = split_routed_mass(forest)))
  if (effects) out$effects <- global_effects(forest, out)
  out
}

# design variants: "joint" (TreeIPF), "marginal" (product of TreeIPF marginals),
# "uniform" (uniform over the padded split box)
design_variant <- function(forest, fit, type = c("joint", "marginal", "uniform")) {
  type <- match.arg(type)
  M <- nrow(fit$Z)
  if (type == "joint") return(list(Z = fit$Z, q = fit$q, mu = fit$mu))
  if (type == "marginal") {
    Zp <- fit$Z
    for (j in seq_len(ncol(Zp))) Zp[, j] <- sample(fit$Z[, j], M, TRUE, prob = fit$q)
  } else {
    Zp <- build_grid(forest, M, fit$grid$sets, design = "uniform")$Z
  }
  list(Z = Zp, q = rep(1 / M, M), mu = forest_mean(forest, Zp))
}

# the 3 x 2 x 2 factorial: design x trimming x SE scaling
ablation_grid <- function(forest, fit, alpha = 0.05) {
  res <- list()
  for (dt in c("uniform", "marginal", "joint")) {
    d <- design_variant(forest, fit, dt)
    for (tm in c(TRUE, FALSE)) for (sm in c("rootn", "grid")) {
      pr <- project_design(forest, d$Z, d$q, d$mu, fit$grid, alpha, trim = tm, se_mode = sm)
      key <- paste(dt, if (tm) "trim" else "full", sm, sep = "|")
      res[[key]] <- pr
    }
  }
  res
}

# partial-dependence and accumulated-local-effects average slopes on the link scale,
# computed over the reconstructed law Q* and restricted to the same interior
global_effects <- function(forest, fit, n_bg = 300L, n_pts = 15L, n_bins = 10L) {
  Z <- fit$Z; q <- fit$q; vt <- fit$grid$var_type; cuts <- fit$cuts
  p <- ncol(Z); sets <- fit$grid$sets
  bg <- Z[sample.int(nrow(Z), n_bg, TRUE, prob = q), , drop = FALSE]
  pd <- ale <- rep(NA_real_, p)
  link <- function(ZZ) forest_link(forest, ZZ)
  for (j in seq_len(p)) {
    if (vt[j] == "unsplit") next
    if (vt[j] == "binary") {
      Z1 <- Z; Z1[, j] <- 1; Z0 <- Z; Z0[, j] <- 0
      ale[j] <- sum(q * (link(Z1) - link(Z0)))
      B1 <- bg; B1[, j] <- 1; B0 <- bg; B0[, j] <- 0
      pd[j] <- mean(link(B1) - link(B0))
      next
    }
    lo <- if (cuts$trimmed[j]) cuts$lo[j] else min(sets[[j]])
    hi <- if (cuts$trimmed[j]) cuts$hi[j] else max(sets[[j]])
    sel <- Z[, j] >= lo & Z[, j] <= hi
    if (sum(sel) < 10) next
    # PD: slope of the PD curve over weighted quantiles of Q*_j
    o <- order(Z[sel, j]); cw <- cumsum(q[sel][o]) / sum(q[sel])
    g <- unique(Z[sel, j][o][pmin(findInterval(seq(0.5, n_pts - 0.5) / n_pts, cw) + 1L, sum(sel))])
    pdv <- vapply(g, function(v) { Bv <- bg; Bv[, j] <- v; mean(link(Bv)) }, 0)
    pd[j] <- if (length(g) > 1) coef(lm(pdv ~ g))[2] else NA_real_
    # ALE: equal-mass bins of Q*_j inside the interior, mass-weighted local slopes
    edges <- unique(c(lo, Z[sel, j][o][pmin(findInterval((1:(n_bins - 1)) / n_bins, cw) + 1L, sum(sel))], hi))
    if (length(edges) < 2) next
    zs <- Z[sel, , drop = FALSE]; qs <- q[sel]
    k <- pmin(findInterval(zs[, j], edges, rightmost.closed = TRUE), length(edges) - 1L)
    Zu <- zs; Zu[, j] <- edges[k + 1L]; Zl <- zs; Zl[, j] <- edges[k]
    slope <- (link(Zu) - link(Zl)) / (edges[k + 1L] - edges[k])
    ale[j] <- sum(qs * slope) / sum(qs)
  }
  list(pd = setNames(c(NA, pd), c("(Intercept)", forest$var_names)),
       ale = setNames(c(NA, ale), c("(Intercept)", forest$var_names)))
}

# ------------------------------------------------------------------------------
# 8. parametric bootstrap that retrains with the released hyperparameters
# ------------------------------------------------------------------------------

tree_ipf_boot <- function(forest, fit, B_boot = 30L, M = 2200L, alpha = 0.05) {
  n <- forest$n_hat
  keep <- c("est", "se")
  draws <- matrix(NA_real_, B_boot, length(fit$est), dimnames = list(NULL, names(fit$est)))
  vn <- forest$var_names
  for (r in seq_len(B_boot)) {
    idx <- sample.int(nrow(fit$Z), n, TRUE, prob = fit$q)
    Xs <- fit$Z[idx, , drop = FALSE]
    m <- fit$mu[idx]
    ys <- if (forest$family == "gaussian") {
      s2 <- if (!is.na(fit$sigma2)) fit$sigma2 else forest$sigma2
      m + rnorm(n, 0, sqrt(s2))
    } else rbinom(n, 1, m)
    fs <- tryCatch(refit_like(forest, Xs, ys), error = function(e) NULL)
    if (is.null(fs)) next
    fb <- tryCatch(tree_ipf(fs, M = M, alpha = alpha), error = function(e) NULL)
    if (!is.null(fb)) draws[r, ] <- fb$est
  }
  list(draws = draws, se = apply(draws, 2, sd, na.rm = TRUE),
       bias = colMeans(draws, na.rm = TRUE) - fit$est)
}

# ------------------------------------------------------------------------------
# 9. reference fits that need the patient data (simulation only)
# ------------------------------------------------------------------------------

ipd_glm <- function(X, y, family, rows = NULL) {
  if (!is.null(rows)) { X <- X[rows, , drop = FALSE]; y <- y[rows] }
  d <- data.frame(y = y, X)
  f <- if (family == "gaussian") lm(y ~ ., d) else glm(y ~ ., d, family = binomial())
  cf <- summary(f)$coefficients
  est <- se <- setNames(rep(NA_real_, ncol(X) + 1L), c("(Intercept)", colnames(X)))
  est[rownames(cf)] <- cf[, 1]; se[rownames(cf)] <- cf[, 2]
  list(est = est, se = se)
}

# the ensemble projected over the cohort's own covariates (perfect step 1)
fhat_on_true_x <- function(forest, X, cuts) {
  rows <- in_box(X, cuts)
  mu <- forest_mean(forest, X)
  vt <- var_types(split_sets(forest)$sets, at_values(forest))
  keep <- which(vt != "unsplit")
  D <- cbind(1, X[rows, keep, drop = FALSE])
  b <- wls_glm(D, rep(1 / sum(rows), sum(rows)), mu[rows], forest$family)
  est <- setNames(rep(NA_real_, ncol(X) + 1L), c("(Intercept)", colnames(X)))
  est[c(1L, keep + 1L)] <- b
  est
}
