# ==============================================================================
# 10_summarise_results.R  --  collect replicate files and compute performance
# measures with Monte Carlo SEs (Morris et al. 2019).  Unconditional summaries:
# a coefficient that could not be estimated (predictor never split) counts as
# not covering and not rejecting.
# ==============================================================================

suppressPackageStartupMessages(library(data.table))
root_dir <- normalizePath(Sys.getenv("TREEIPF_ROOT", "."), winslash = "/")
res_dir <- file.path(root_dir, "results")

collect <- function(study, part = "est") {
  fs <- list.files(file.path(res_dir, "reps", study), full.names = TRUE, pattern = "\\.rds$")
  if (!length(fs)) return(NULL)
  r <- lapply(fs, readRDS)
  bad <- vapply(r, function(x) is.null(x) || !is.null(x$error), TRUE)
  if (any(bad)) message(study, ": ", sum(bad), " replicate files with errors")
  r <- r[!bad]
  if (is.data.table(r[[1]])) rbindlist(r, fill = TRUE) else rbindlist(lapply(r, `[[`, part), fill = TRUE)
}

perf <- function(d, by) {
  d[, {
    ok <- is.finite(est)
    R <- .N; e <- est[ok] - target[ok]
    has_se <- is.finite(se)
    lo <- est - 1.96 * se; hi <- est + 1.96 * se
    cov <- ifelse(has_se & ok, lo <= target & target <= hi, FALSE)
    rej <- ifelse(has_se & ok, abs(est / se) > 1.96, FALSE)
    sdv <- sd(est[ok])
    list(R = R, n_na = sum(!ok),
         mean_est = mean(est[ok]), bias = mean(e), mcse_bias = sd(e) / sqrt(sum(ok)),
         rel_bias = mean(e) / ifelse(target[1] == 0, NA_real_, target[1]),
         emp_sd = sdv, mean_se = mean(se[has_se & ok]),
         se_sd = mean(se[has_se & ok]) / sdv,
         mcse_se_sd = (mean(se[has_se & ok]) / sdv) / sqrt(2 * (sum(ok) - 1)),
         rmse = sqrt(mean(e^2)),
         coverage = if (any(has_se)) mean(cov) else NA_real_,
         mcse_cov = if (any(has_se)) sqrt(mean(cov) * (1 - mean(cov)) / R) else NA_real_,
         reject = if (any(has_se)) mean(rej) else NA_real_,
         mcse_rej = if (any(has_se)) sqrt(mean(rej) * (1 - mean(rej)) / R) else NA_real_,
         target = mean(target), target_sd = sd(target))
  }, by = by]
}

keys <- c("study", "family", "scenario", "n", "rho", "depth", "alpha", "b5", "method", "term")

summarise_study <- function(study) {
  d <- collect(study)
  if (is.null(d)) return(NULL)
  s <- perf(d, keys)
  # bias decomposition: ensemble part and reconstruction part (eq. 4)
  w <- dcast(d[method %in% c("TreeIPF", "Ensemble on true X")],
             study + family + scenario + n + rho + depth + alpha + b5 + rep + term ~ method,
             value.var = c("est", "target"))
  setnames(w, c("est_Ensemble on true X", "est_TreeIPF", "target_TreeIPF"), c("ens", "ipf", "tgt"), skip_absent = TRUE)
  if (all(c("ens", "ipf", "tgt") %in% names(w))) {
    dec <- w[, .(ensemble_part = mean(ens - tgt, na.rm = TRUE),
                 reconstruction_part = mean(ipf - ens, na.rm = TRUE),
                 sd_reconstruction = sd(ipf - ens, na.rm = TRUE)),
             by = .(study, family, scenario, n, rho, depth, alpha, b5, term)]
  } else dec <- NULL
  info <- if ("rho12" %in% names(d)) d[method == "TreeIPF", .(rho12 = mean(rho12, na.rm = TRUE),
          w_alpha = mean(w_alpha, na.rm = TRUE), kl = mean(kl, na.rm = TRUE), iter = mean(iter, na.rm = TRUE),
          B = mean(B, na.rm = TRUE), q_binds = mean(q_binds > 0, na.rm = TRUE), fallback = mean(fallback, na.rm = TRUE)),
          by = .(study, family, scenario, n, rho, depth, alpha, b5)] else NULL
  list(perf = s, decomp = dec, info = info)
}

if (sys.nframe() == 0L) {
  for (st in c("S1", "S2", "S3", "S4", "S5")) {
    s <- summarise_study(st)
    if (is.null(s)) next
    fwrite(s$perf, file.path(res_dir, sprintf("sim_%s_performance.csv", st)))
    if (!is.null(s$decomp)) fwrite(s$decomp, file.path(res_dir, sprintf("sim_%s_decomposition.csv", st)))
    if (!is.null(s$info)) fwrite(s$info, file.path(res_dir, sprintf("sim_%s_info.csv", st)))
    dg <- collect(st, "diag")
    if (!is.null(dg)) fwrite(dg[, .(n_splits = mean(n_splits), unsplit = mean(n_splits == 0),
                                    routed_mass = mean(routed_mass)),
                                by = .(study, family, scenario, n, rho, depth, b5, term)],
                             file.path(res_dir, sprintf("sim_%s_splits.csv", st)))
    cat(st, "summarised\n")
  }
  # Study 6: bootstrap
  b <- collect("S6")
  if (!is.null(b)) {
    b <- b[term != "(Intercept)"]
    b[, se_use := pmax(se_boot, se_model)]
    b[, est_bc := est - bias_boot]
    sm <- b[, {
      cv <- function(lo, hi) { x <- lo <= target & target <= hi; x[is.na(x)] <- FALSE; c(mean(x), sqrt(mean(x) * (1 - mean(x)) / .N)) }
      a <- cv(est - 1.96 * se_model, est + 1.96 * se_model)
      bb <- cv(est - 1.96 * se_boot, est + 1.96 * se_boot)
      bc <- cv(est_bc - 1.96 * se_use, est_bc + 1.96 * se_use)
      ba <- cv(2 * est - q975, 2 * est - q025)
      list(R = .N, bias = mean(est - target), bias_bc = mean(est_bc - target, na.rm = TRUE),
           emp_sd = sd(est), emp_sd_bc = sd(est_bc, na.rm = TRUE),
           mean_se_model = mean(se_model), mean_se_boot = mean(se_boot, na.rm = TRUE),
           rmse = sqrt(mean((est - target)^2)), rmse_bc = sqrt(mean((est_bc - target)^2, na.rm = TRUE)),
           cov_model = a[1], cov_boot = bb[1], cov_bc = bc[1], cov_basic = ba[1], mcse_cov_bc = bc[2],
           n_ok = mean(n_ok))
    }, by = .(family, scenario, n, rho, term)]
    fwrite(sm, file.path(res_dir, "sim_S6_bootstrap.csv"))
    cat("S6 summarised\n")
  }
  # clinical audits: Lin's concordance on standardised coefficients (1/n moments),
  # agreement per run, and 50-run medians of the rebuilt moments
  cf <- if (file.exists(file.path(res_dir, "clinical_audit_coefficients.csv"))) fread(file.path(res_dir, "clinical_audit_coefficients.csv")) else NULL
  if (!is.null(cf)) {
    mo <- fread(file.path(res_dir, "clinical_audit_moments.csv")); cr <- fread(file.path(res_dir, "clinical_audit_correlations.csv"))
    lin <- function(x, y) { mx <- mean(x); my <- mean(y)
      2 * mean((x - mx) * (y - my)) / (mean((x - mx)^2) + mean((y - my)^2) + (mx - my)^2) }
    sdv <- mo[seed == 1, .(cohort, term, emp_sd, bin = var_type == "binary")]
    c2 <- merge(cf[term != "(Intercept)"], sdv, by = c("cohort", "term"))
    c2[, std := ifelse(bin, est, est * emp_sd)]
    g <- c2[method == "glm", .(cohort, seed, term, gstd = std, gest = est, gse = se, gz = est / se)]
    d <- merge(c2[method != "glm"], g, by = c("cohort", "seed", "term"))
    per_run <- d[, .(ccc = lin(gstd, std), inside = sum(abs(est - gest) <= 1.96 * se, na.rm = TRUE), k = .N,
                     se_ratio = median(se / gse, na.rm = TRUE),
                     lost_sig = sum(abs(gz) > 1.96 & abs(est / se) <= 1.96, na.rm = TRUE), glm_sig = sum(abs(gz) > 1.96)),
                 by = .(cohort, seed, method)]
    fwrite(per_run[, .(ccc_median = median(ccc), ccc_p05 = quantile(ccc, 0.05), ccc_p95 = quantile(ccc, 0.95),
                       inside_median = median(inside), k = k[1], se_ratio_median = median(se_ratio),
                       lost_sig_median = median(lost_sig), glm_sig = glm_sig[1]), by = .(cohort, method)],
           file.path(res_dir, "clinical_audit_summary.csv"))
    fwrite(d[, .(median_est = median(est), p05 = quantile(est, 0.05), p95 = quantile(est, 0.95), glm = gest[1],
                 share_negative = mean(est < 0), share_sig_neg = mean(est / se < -1.96, na.rm = TRUE),
                 share_sig_pos = mean(est / se > 1.96, na.rm = TRUE), share_inside = mean(abs(est - gest) <= 1.96 * se, na.rm = TRUE)),
             by = .(cohort, method, term)], file.path(res_dir, "clinical_audit_terms.csv"))
    fwrite(mo[, .(emp_mean = emp_mean[1], rec_mean = median(rec_mean), emp_sd = emp_sd[1], rec_sd = median(rec_sd),
                  n_splits = median(n_splits), routed_mass = median(routed_mass)), by = .(cohort, term)],
           file.path(res_dir, "clinical_moments_summary.csv"))
    fwrite(cr[, .(r = cor(emp, rec), mae = mean(abs(emp - rec))), by = .(cohort, seed)][
             , .(r_median = median(r), r_min = min(r), r_max = max(r), mae_median = median(mae)), by = cohort],
           file.path(res_dir, "clinical_correlation_summary.csv"))
    cat("clinical summarised\n")
  }
  l <- collect("S7")
  if (!is.null(l)) {
    fwrite(perf(l, c("method", "family", "term")), file.path(res_dir, "sim_S7_libraries.csv"))
    fwrite(l[!is.na(n_hat), .(n_hat_mean = mean(n_hat), n_hat_sd = sd(n_hat), sigma2 = mean(sigma2, na.rm = TRUE),
                              rho12 = mean(rho12), miss_mass = mean(miss_mass), B = mean(B)),
             by = .(method, family)], file.path(res_dir, "sim_S7_library_info.csv"))
    cat("S7 summarised\n")
  }
}
