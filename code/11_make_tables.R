# ==============================================================================
# 11_make_tables.R  --  LaTeX tables (manuscript and Additional file 1) built
# directly from the result CSVs, so that no number is typed by hand.
# ==============================================================================

suppressPackageStartupMessages(library(data.table))
root_dir <- normalizePath(Sys.getenv("TREEIPF_ROOT", "."), winslash = "/")
res <- file.path(root_dir, "results")
tab_dir <- file.path(root_dir, "manuscript", "tables")
dir.create(tab_dir, showWarnings = FALSE)

f3 <- function(x, d = 3) ifelse(is.na(x), "--", formatC(x, format = "f", digits = d))
pct <- function(x, d = 0) ifelse(is.na(x), "--", paste0(formatC(100 * x, format = "f", digits = d)))
neg <- function(s) gsub("^-", "$-$", s)
rd <- function(f) { p <- file.path(res, f); if (file.exists(p)) fread(p) else NULL }
write_tab <- function(lines, name) writeLines(lines, file.path(tab_dir, name))
fam_lab <- c(gaussian = "Linear", binomial = "Logistic")

# ---- Table 2: TreeIPF across Study 1 -------------------------------------------
tab_study1 <- function() {
  p <- rd("sim_S1_performance.csv"); inf <- rd("sim_S1_info.csv")
  if (is.null(p)) return(invisible())
  g <- p[method == "TreeIPF"]
  w <- dcast(g, family + n + rho ~ term, value.var = c("rel_bias", "se_sd", "coverage", "reject"))
  w <- merge(w, inf[, .(family, n, rho, rho12)], by = c("family", "n", "rho"))
  w[, family := factor(family, c("gaussian", "binomial"))]; setorder(w, family, rho, n)
  L <- c("\\begin{table*}[!t]",
         "\\caption{\\TreeIPF{} in Study 1: relative bias, SE/SD ratio and coverage for three predictors, and Type I error for the null predictor}\\label{tab:study1}",
         "\\footnotesize\\setlength{\\tabcolsep}{2.6pt}",
         "\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}llrrrrrrrrrrr@{}}", "\\toprule",
         " & & & \\multicolumn{3}{c}{$X_1$ (correlated)} & \\multicolumn{3}{c}{$X_3$ (skewed)} & \\multicolumn{3}{c}{$X_4$ (binary)} & $X_5$ (null) \\\\",
         "\\cmidrule(lr){4-6}\\cmidrule(lr){7-9}\\cmidrule(lr){10-12}\\cmidrule(l){13-13}",
         "Model & $\\rho$ & $n$ & Bias \\% & SE/SD & Cover & Bias \\% & SE/SD & Cover & Bias \\% & SE/SD & Cover & Type I \\\\", "\\midrule")
  for (i in seq_len(nrow(w))) {
    r <- w[i]
    L <- c(L, paste0(if (i == 1 || r$family != w$family[i - 1]) fam_lab[as.character(r$family)] else "", " & ",
                     if (i == 1 || r$rho != w$rho[i - 1] || r$family != w$family[i - 1]) f3(r$rho, 1) else "", " & ", r$n, " & ",
                     neg(pct(r$rel_bias_X1)), " & ", f3(r$se_sd_X1, 2), " & ", f3(r$coverage_X1, 2), " & ",
                     neg(pct(r$rel_bias_X3)), " & ", f3(r$se_sd_X3, 2), " & ", f3(r$coverage_X3, 2), " & ",
                     neg(pct(r$rel_bias_X4)), " & ", f3(r$se_sd_X4, 2), " & ", f3(r$coverage_X4, 2), " & ",
                     f3(r$reject_X5, 3), " \\\\"))
    if (i < nrow(w) && w$family[i + 1] != r$family) L <- c(L, "\\midrule")
  }
  L <- c(L, "\\botrule", "\\end{tabular*}",
         "\\par\\vspace{2pt}\\noindent{\\scriptsize Bias \\%: mean relative bias against the target (here the true coefficient). SE/SD: mean model-based SE divided by the empirical SD of the estimates. Cover: coverage of nominal 95\\% intervals. 200 replicates per row. Largest Monte Carlo SE in the table: ", paste0(f3(max(g$mcse_cov[g$term != "(Intercept)"], na.rm = TRUE), 3), " for coverage, ", pct(max(g[term %in% c("X1","X3","X4"), mcse_bias / abs(target)]), 1), "\\% for relative bias, ", f3(max(g$mcse_rej[g$term == "X5"]), 3), " for Type I error. Other methods and the rebuilt correlations: Additional file~1, Tables~S1--S3.}"),
         "\\end{table*}")
  write_tab(L, "tab_study1.tex")
}

# ---- Table 3: library files (Study 7) ---------------------------------------------
tab_libs <- function() {
  p <- rd("sim_S7_libraries.csv"); inf <- rd("sim_S7_library_info.csv")
  if (is.null(p)) return(invisible())
  m <- c("Oracle", "TreeIPF XGBoost", "TreeIPF XGBoost (covers as counts)", "TreeIPF sklearn", "TreeIPF ranger",
         "Unif-rootn XGBoost", "Unif-rootn sklearn", "Unif-rootn ranger")
  lab <- c("Oracle (full data)", "XGBoost: TreeIPF", "XGBoost: covers used as counts", "scikit-learn: TreeIPF",
           "ranger: TreeIPF", "XGBoost: Unif-rootn", "scikit-learn: Unif-rootn", "ranger: Unif-rootn")
  L <- c("\\begin{table*}[!t]",
         "\\caption{\\TreeIPF{} on model files saved by XGBoost, scikit-learn and ranger (Study 7; $n=1500$, $\\rho=0.4$, 100 cohorts)}\\label{tab:libs}",
         "\\footnotesize\\setlength{\\tabcolsep}{3.5pt}",
         "\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}llrrrrrrr@{}}", "\\toprule",
         " & & & \\multicolumn{2}{c}{$X_1$ (correlated)} & \\multicolumn{2}{c}{$X_4$ (binary)} & $X_5$ & \\\\",
         "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-8}",
         "Model & File and method & $\\hat n$ & Bias \\% & Cover & Bias \\% & Cover & Type I & SE/SD ($X_1$) \\\\", "\\midrule")
  for (fam in c("gaussian", "binomial")) {
    for (k in seq_along(m)) {
      e <- p[family == fam & method == m[k]]
      if (!nrow(e)) next
      x1 <- e[term == "X1"]; x4 <- e[term == "X4"]; x5 <- e[term == "X5"]
      nh <- inf[family == fam & method == m[k]]$n_hat_mean
      L <- c(L, paste0(if (k == 1) fam_lab[fam] else "", " & ", lab[k], " & ", if (length(nh)) formatC(nh, format = "f", digits = 0) else "--", " & ",
                       neg(pct(x1$rel_bias)), " & ", f3(x1$coverage, 2), " & ", neg(pct(x4$rel_bias)), " & ", f3(x4$coverage, 2), " & ",
                       f3(x5$reject, 3), " & ", f3(x1$se_sd, 2), " \\\\"))
    }
    if (fam == "gaussian") L <- c(L, "\\midrule")
  }
  L <- c(L, "\\botrule", "\\end{tabular*}",
         "\\par\\vspace{2pt}\\noindent{\\scriptsize $\\hat n$: cohort size read from the file (true $n=1500$). Bias \\%, Cover, Type I and SE/SD as in Table~\\ref{tab:study1}; Monte Carlo SE about 0.022 for coverage near 0.95 with 100 cohorts.}",
         "\\end{table*}")
  write_tab(L, "tab_libs.tex")
}

# ---- Table 4: IST pooled odds ratios ------------------------------------------------
tab_ist <- function() {
  ma <- rd("ist_meta_analysis.csv")
  if (is.null(ma)) return(invisible())
  lab <- c(age10 = "Age (per 10 years)", female = "Female sex", sbp10 = "Systolic BP (per 10 mmHg)",
           consc_drowsy = "Impaired consciousness", ratrial = "Atrial fibrillation", ct_inf = "Infarct on CT",
           rdef1 = "Facial deficit", rdef2 = "Arm or hand deficit")
  st <- c("Gold: GLM in all 12", "All 12 from ensembles: TreeIPF", "All 12 from ensembles: TreeIPF + bootstrap",
          "All 12 from ensembles: Unif-rootn", "All 12 from ensembles: Naive-Unif")
  cell <- function(e) sprintf("%s (%s, %s)", f3(e$OR, 2), f3(e$OR_lb, 2), f3(e$OR_ub, 2))
  L <- c("\\begin{table*}[!t]",
         "\\caption{International Stroke Trial: pooled odds ratios (95\\% CI) and $I^2$ when all 12 countries enter from their full-data regressions or from their boosted ensembles}\\label{tab:ist}",
         "\\footnotesize\\setlength{\\tabcolsep}{3pt}",
         "\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}lrrrrrr@{}}", "\\toprule",
         " & \\multicolumn{2}{c}{Full-data regressions} & \\multicolumn{2}{c}{\\TreeIPF{}} & \\multicolumn{2}{c}{\\TreeIPF{} + bootstrap} \\\\",
         "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(l){6-7}",
         "Predictor & OR (95\\% CI) & $I^2$ & OR (95\\% CI) & $I^2$ & OR (95\\% CI) & $I^2$ \\\\", "\\midrule")
  for (tm in names(lab)) {
    g <- ma[term == tm & strategy == st[1]]; t <- ma[term == tm & strategy == st[2]]
    b <- ma[term == tm & strategy == st[3]]
    L <- c(L, paste0(lab[tm], " & ", cell(g), " & ", f3(g$I2, 0), " & ", cell(t), " & ", f3(t$I2, 0), " & ",
                     cell(b), " & ", f3(b$I2, 0), " \\\\"))
  }
  L <- c(L, "\\botrule", "\\end{tabular*}",
         "\\par\\vspace{2pt}\\noindent{\\scriptsize Random-effects meta-analysis (REML) of 12 countries (15{,}580 patients). Full-data regressions: logistic regression on each country's data. Other columns: each country enters only through its boosted ensemble. $I^2$ in \\%. Naive-Unif and Unif-rootn: Fig.~\\ref{fig:ist}; hybrid analyses: Additional file~1.}",
         "\\end{table*}")
  write_tab(L, "tab_ist.tex")
}

# ---- Additional file tables ----------------------------------------------------------
longtab <- function(d, cols, heads, caption, label, name, align) {
  d <- copy(d)
  for (cc in cols) if (is.character(d[[cc]])) d[[cc]] <- ifelse(grepl("^\\$", d[[cc]]), d[[cc]], gsub("_", "\\_", d[[cc]], fixed = TRUE))
  L <- c(sprintf("\\begin{longtable}{%s}", align), sprintf("\\caption{%s}\\label{%s}\\\\", caption, label), "\\toprule",
         paste(heads, collapse = " & "), "\\\\ \\midrule \\endfirsthead", "\\toprule", paste(heads, collapse = " & "),
         "\\\\ \\midrule \\endhead", "\\bottomrule \\endfoot")
  for (i in seq_len(nrow(d))) L <- c(L, paste0(paste(vapply(cols, function(cc) as.character(d[[cc]][i]), ""), collapse = " & "), " \\\\"))
  write_tab(c(L, "\\end{longtable}"), name)
}

tab_s1 <- function() {
  p <- rd("sim_S1_performance.csv"); if (is.null(p)) return(invisible())
  mm <- c(Oracle = "Oracle", `Oracle-interior` = "Oracle inside", `Ensemble on true X` = "Ensemble, real X", TreeIPF = "TreeIPF",
          `TreeIPF-Bag` = "TreeIPF-Bag", `A:uniform|full|grid` = "Naive-Unif", `A:uniform|trim|rootn` = "Unif-rootn",
          `A:marginal|trim|rootn` = "Marg-IPF", ALE = "ALE", PD = "PD")
  d <- p[method %in% names(mm) & term != "(Intercept)"]
  d[, meth := mm[method]]; d[, family := fam_lab[family]]
  setorder(d, family, n, rho, term, method)
  d[, `:=`(b = neg(f3(bias)), s = f3(se_sd, 2), c = f3(coverage, 3), r = f3(reject, 3), mse = f3(rmse))]
  longtab(d, c("family", "n", "rho", "term", "meth", "b", "mse", "s", "c", "r"),
          c("Model", "$n$", "$\\rho$", "Term", "Method", "Bias", "RMSE", "SE/SD", "Cover", "Reject"),
          "Study 1, all settings and methods. Reject: rejection rate of the 5\\% Wald test (Type I error for $X_5$, power otherwise).",
          "tab:s1", "tab_s1_study1.tex", "llrllrrrrr")
}

tab_s2 <- function() {
  p <- rd("sim_S5_performance.csv"); if (is.null(p)) return(invisible())
  mm <- c(Oracle = "Oracle (full-range target)", `Oracle-interior` = "Oracle inside", TreeIPF = "TreeIPF",
          `TreeIPF vs full target` = "TreeIPF vs full-range target", `A:uniform|trim|rootn` = "Unif-rootn",
          `A:marginal|trim|rootn` = "Marg-IPF", PD = "PD")
  d <- p[method %in% names(mm) & term %in% c("X1", "X2", "X3", "X4", "X5")]
  d[, meth := mm[method]]; d[, family := fam_lab[family]]
  setorder(d, scenario, family, n, term, method)
  d[, `:=`(t = f3(target), b = neg(f3(bias)), s = f3(se_sd, 2), c = f3(coverage, 3), r = f3(reject, 3), na = n_na)]
  longtab(d, c("scenario", "family", "n", "term", "meth", "t", "b", "s", "c", "r", "na"),
          c("Setting", "Model", "$n$", "Term", "Method", "Target", "Bias", "SE/SD", "Cover", "Reject", "Not split"),
          "Study 5: wrong working model (nonlinear), rare events and $p=12$. Target: $\\theta_{\\mathcal I}$ unless stated.",
          "tab:s2", "tab_s2_study5.tex", "lllllrrrrrr")
}

tab_s3 <- function() {
  b <- rd("sim_S6_bootstrap.csv"); if (is.null(b)) return(invisible())
  b[, family := fam_lab[family]]; setorder(b, scenario, family, n, rho, term)
  b[, `:=`(b0 = neg(f3(bias)), b1 = neg(f3(bias_bc)), c0 = f3(cov_model, 2), c1 = f3(cov_boot, 2), c2 = f3(cov_bc, 2), c3 = f3(cov_basic, 2),
           s0 = f3(mean_se_model), s1 = f3(mean_se_boot), sd0 = f3(emp_sd))]
  longtab(b, c("scenario", "family", "n", "rho", "term", "b0", "b1", "sd0", "s0", "s1", "c0", "c1", "c2", "c3"),
          c("Setting", "Model", "$n$", "$\\rho$", "Term", "Bias", "Bias (corr.)", "Emp SD", "SE model", "SE boot", "Cov model", "Cov boot", "Cov corr.", "Cov basic"),
          "Study 6: parametric bootstrap (30 retrained ensembles per cohort, 100 cohorts per setting). Cov corr.: bias-corrected estimate with the larger of the two SEs; Cov basic: basic bootstrap interval.",
          "tab:s3", "tab_s3_bootstrap.tex", "lllllrrrrrrrrr")
}

tab_s4 <- function() {
  p2 <- rd("sim_S2_performance.csv"); p3 <- rd("sim_S3_performance.csv"); i2 <- rd("sim_S2_info.csv")
  if (is.null(p2) || is.null(p3)) return(invisible())
  d2 <- p2[method %in% c("TreeIPF", "Marg-IPF") & term %in% c("X1", "X4", "X5")]
  d2 <- merge(d2, i2[, .(family, depth, rho12)], by = c("family", "depth"))
  d2[, family := fam_lab[family]]; setorder(d2, family, depth, term, method)
  d2[, `:=`(b = neg(f3(bias)), s = f3(se_sd, 2), c = f3(coverage, 3), r = f3(reject, 3), rr = f3(rho12, 2), na = n_na)]
  longtab(d2, c("family", "depth", "term", "method", "rr", "b", "s", "c", "r", "na"),
          c("Model", "Depth", "Term", "Method", "Rebuilt $\\rho$", "Bias", "SE/SD", "Cover", "Reject", "Not split"),
          "Study 2: tree depth ($n=800$, $\\rho=0.7$). Not split: replicates (of 200) in which the predictor was never split.",
          "tab:s4a", "tab_s4_depth_trim.tex", "llllrrrrrr")
  d3 <- p3[grepl("alpha=", method) & term %in% c("X1", "X3", "X4")]
  d3[, a := sub("TreeIPF alpha=", "", method)]; d3[, family := fam_lab[family]]; setorder(d3, family, term, a)
  d3[, `:=`(b = neg(f3(bias)), s = f3(se_sd, 2), c = f3(coverage, 3), m = f3(mean_se), sdv = f3(emp_sd))]
  tmp <- file.path(tab_dir, "tmp.tex")
  longtab(d3, c("family", "term", "a", "b", "sdv", "m", "s", "c"),
          c("Model", "Term", "$\\alpha$", "Bias", "Emp SD", "Mean SE", "SE/SD", "Cover"),
          "Study 3: trimming budget $\\alpha$ applied to the same fitted ensembles ($n=800$, $\\rho=0.4$).",
          "tab:s4b", "tmp.tex", "lllrrrrr")
  write_tab(c(readLines(file.path(tab_dir, "tab_s4_depth_trim.tex")), readLines(tmp)), "tab_s4_depth_trim.tex"); file.remove(tmp)
}

tab_s5 <- function() {
  cf <- rd("clinical_audit_coefficients.csv"); if (is.null(cf)) return(invisible())
  g <- cf[method == "glm", .(glm = mean(est), glm_se = mean(se)), by = .(cohort, term)]
  d <- merge(cf[method != "glm" & term != "(Intercept)"], g, by = c("cohort", "term"))
  s <- d[, .(diff = median(est - glm), se_ratio = median(se / glm_se, na.rm = TRUE),
             inside = mean(abs(est - glm) <= 1.96 * se, na.rm = TRUE)), by = .(cohort, method, term)]
  s <- s[method %in% c("TreeIPF", "Unif-rootn", "Naive-Unif", "Marg-IPF")]
  setorder(s, cohort, term, method)
  s[, `:=`(d = neg(f3(diff)), r = f3(se_ratio, 2), i = f3(inside, 2))]
  longtab(s, c("cohort", "term", "method", "d", "r", "i"),
          c("Data", "Predictor", "Method", "Median difference", "Median SE ratio", "Share inside"),
          "Clinical audits over 50 seeds: median difference from the full-data logistic coefficient, median ratio of SEs, and share of seeds in which the full-data estimate lies inside the method's 95\\% interval.",
          "tab:s5", "tab_s5_clinical.tex", "lllrrr")
}

tab_s6 <- function() {
  ce <- rd("ist_country_estimates.csv"); if (is.null(ce)) return(invisible())
  d <- ce[term != "(Intercept)"]
  d[, `:=`(diff = ipf_est - glm_est, ratio = ipf_se / glm_se, inside = abs(ipf_est - glm_est) <= 1.96 * ipf_se)]
  s <- d[, .(n = n[1], mean_diff = mean(diff), med_ratio = median(ratio), inside = sum(inside), shrunk = sum(abs(ipf_est) < abs(glm_est))),
         by = .(rank_id, country_code)][order(rank_id)]
  s[, `:=`(md = neg(f3(mean_diff)), mr = f3(med_ratio, 2), ins = paste0(inside, "/8"), shr = paste0(shrunk, "/8"))]
  longtab(s, c("country_code", "n", "md", "mr", "ins", "shr"),
          c("Country", "$n$", "Mean difference", "Median SE ratio", "Inside CI", "Closer to 0"),
          "International Stroke Trial, country level: TreeIPF from each country's ensemble against the same country's full-data regression (8 predictors).",
          "tab:s6", "tab_s6_ist_country.tex", "lrrrrr")
}


tab_s7 <- function() {
  ma <- rd("ist_meta_analysis.csv"); mr <- rd("ist_meta_regression.csv"); if (is.null(ma)) return(invisible())
  d <- ma[grepl("Hybrid|Classical", strategy)]
  d[, `:=`(o = sprintf("%s (%s, %s)", f3(OR, 2), f3(OR_lb, 2), f3(OR_ub, 2)), i = f3(I2, 0), w = f3(ml_weight, 0))]
  setorder(d, term, strategy)
  longtab(d, c("term", "strategy", "o", "i", "w"), c("Predictor", "Analysis", "OR (95\\% CI)", "$I^2$", "Weight of ML countries (\\%)"),
          "International Stroke Trial: hybrid analyses (six countries from regressions, six from ensembles, both splits) and analyses of six regression countries only.",
          "tab:s7", "tab_s7_hybrid.tex", "lllrr")
  if (!is.null(mr)) {
    mr[, `:=`(sh = neg(f3(ml_shift)), se = f3(ml_shift_se), pp = formatC(p, format = "g", digits = 2))]
    L <- readLines(file.path(tab_dir, "tab_s7_hybrid.tex"))
    tmp <- "tmp2.tex"
    longtab(mr, c("term", "sh", "se", "pp"), c("Predictor", "Shift (log OR)", "SE", "$p$"),
            "Meta-regression in the hybrid analysis (odd-ranked countries from ensembles): shift of ensemble-based estimates relative to regression-based ones.",
            "tab:s7b", tmp, "lrrr")
    write_tab(c(L, readLines(file.path(tab_dir, tmp))), "tab_s7_hybrid.tex"); file.remove(file.path(tab_dir, tmp))
  }
}

tab_s1b <- function() {
  p <- rd("sim_S1_performance.csv"); inf <- rd("sim_S1_info.csv"); if (is.null(p)) return(invisible())
  a <- p[grepl("^A:", method)]
  a[, c("design", "region", "se") := tstrsplit(sub("A:", "", method), "|", fixed = TRUE)]
  s <- merge(a[term %in% c("X1", "X2", "X3", "X4"), .(rb = median(rel_bias), cov = median(coverage), sesd = median(se_sd), rmse = median(rmse)),
               by = .(family, design, region, se)],
             a[term == "X5", .(t1lo = min(reject), t1hi = max(reject)), by = .(family, design, region, se)],
             by = c("family", "design", "region", "se"))
  s[, family := fam_lab[family]]
  s[, design := c(uniform = "uniform", marginal = "product of marginals", joint = "joint (TreeIPF)")[design]]
  s[, region := c(trim = "inside outer splits", full = "full grid")[region]]
  s[, se := c(rootn = "cohort size", grid = "grid size")[se]]
  setorder(s, family, design, region, se)
  s[, `:=`(b = neg(pct(rb)), c = f3(cov, 2), r = f3(sesd, 2), m = f3(rmse), t = paste0(f3(t1lo, 3), "--", f3(t1hi, 3)))]
  longtab(s, c("family", "design", "region", "se", "b", "m", "r", "c", "t"),
          c("Model", "Distribution", "Region", "SE scaled by", "Bias \\%", "RMSE", "SE/SD", "Cover", "Type I ($X_5$)"),
          "Study 1 ablation: the 12 combinations of predictor distribution, region and SE scaling, applied to the same boosted files. Bias \\%, RMSE, SE/SD and coverage are medians over the 18 settings and the four non-null coefficients; Type I error is the range over the 18 settings.",
          "tab:s1b", "tab_s1b_ablation.tex", "lllllrrrr")
  inf[, family := fam_lab[family]]; setorder(inf, family, rho, n)
  inf[, `:=`(r = f3(rho12, 3), w = f3(w_alpha, 2), k = signif(kl, 2), b = f3(B, 0), qb = f3(q_binds, 2))]
  longtab(inf, c("family", "rho", "n", "r", "b", "w", "k", "qb"),
          c("Model", "$\\rho$", "$n$", "Rebuilt $\\rho$", "Trees", "Mass inside", "Final divergence", "Quantile cut used"),
          "Study 1, information per setting (means over 200 replicates): rebuilt correlation of $X_1$ and $X_2$, number of trees chosen by cross-validation, share of the rebuilt mass inside the region, final divergence of the SMART fit, and share of replicates in which a quantile bound (rather than an outer split) set at least one cut.",
          "tab:s1c", "tab_s1c_info.tex", "lrrrrrrr")
}

tab_checks <- function() {
  g <- rd("checks_grid.csv"); s2 <- rd("checks_sigma2.csv"); ph <- rd("checks_pad_physionet.csv"); it <- rd("checks_iterations.csv")
  if (is.null(g)) return(invisible())
  g[, family := fam_lab[family]]; setorder(g, family, rho, term)
  g[, `:=`(a = f3(grid_sd), b = f3(mean_se), c = f3(ratio, 2))]
  longtab(g, c("family", "rho", "term", "a", "b", "c"), c("Model", "$\\rho$", "Term", "SD over 10 grids", "Mean SE", "Ratio"),
          "Grid Monte Carlo error: SD of the estimate over ten grids for the same file, against the model-based SE (50 cohorts per setting, $n=1500$).",
          "tab:checks1", "tab_checks.tex", "llrrrr")
  L <- readLines(file.path(tab_dir, "tab_checks.tex"))
  s2[, `:=`(a = f3(sigma2_ratio, 3), b = f3(sigma2_alt_ratio, 3), c = f3(rec_sd_x1, 3), d = f3(rec_sd_x3, 3))]
  longtab(s2, c("rho", "pad", "a", "b", "c", "d"), c("$\\rho$", "Pad", "$\\hat\\sigma^2/\\sigma^2$ (root TSS)", "Alternative", "Rebuilt SD $X_1$", "Rebuilt SD $X_3$"),
          "Linear outcome: residual variance relative to the truth (2.25), for the estimator used in the paper and for an alternative that subtracts the variance of the ensemble's predictions; and rebuilt SDs (true 1) for three pad widths. The slopes were identical for all pad widths.",
          "tab:checks2", "tmp3.tex", "rrrrrr")
  L <- c(L, readLines(file.path(tab_dir, "tmp3.tex")))
  ph <- dcast(ph[term != "(Intercept)"], term + emp_sd ~ pad, value.var = "rec_sd")
  setnames(ph, c("term", "emp_sd", "p10", "p25", "p50"))
  ph[, `:=`(e = f3(emp_sd, 2), a = f3(p10, 2), b = f3(p25, 2), c = f3(p50, 2))]
  longtab(ph, c("term", "e", "a", "b", "c"), c("Predictor", "True SD", "Pad 0.10", "Pad 0.25", "Pad 0.50"),
          "PhysioNet: rebuilt SDs for three pad widths (the slopes were identical).", "tab:checks3", "tmp3.tex", "lrrrr")
  L <- c(L, readLines(file.path(tab_dir, "tmp3.tex"))); file.remove(file.path(tab_dir, "tmp3.tex"))
  write_tab(L, "tab_checks.tex")
}

tab_s8 <- function() {
  t1 <- rd("sim_S8_tuning.csv"); t2 <- rd("sim_S8_wrong_hyper.csv"); if (is.null(t1)) return(invisible())
  t1 <- t1[term != "X5"]; t1[, family := fam_lab[family]]; setorder(t1, family, arm, term)
  t1[, `:=`(a = neg(pct(rel_bias, 1)), b = neg(pct(ens_rel, 1)), c = f3(coverage, 2), d = f3(B, 0))]
  longtab(t1, c("family", "arm", "term", "d", "b", "a", "c"),
          c("Model", "Tuning", "Term", "Trees", "Ensemble part (\\%)", "TreeIPF bias (\\%)", "Cover"),
          "Study 8: two other tunings of the shared ensemble (100 cohorts, $n=1500$, $\\rho=0.4$). Fast: learning rate 0.3, 100 trees, depth 6, minimum leaf 5. Slow: learning rate 0.02, up to 400 trees by cross-validation.",
          "tab:s8a", "tab_s8_hyper.tex", "lllrrrr")
  L <- readLines(file.path(tab_dir, "tab_s8_hyper.tex"))
  t2 <- t2[term != "X5"]
  t2[, arm := c(correct = "correct", tau1 = "subsampling assumed 1", lambda0 = "leaf penalty assumed 0", halfB = "half the trees")[arm]]
  t2[, `:=`(a = neg(pct(rel_bias, 1)), b = neg(pct(rel_bias_bc, 1)), c = f3(cov_model, 2), d = f3(cov_bc, 2))]
  longtab(t2, c("arm", "term", "a", "b", "c", "d"),
          c("Hyperparameters used", "Term", "Bias (\\%)", "Corrected bias (\\%)", "Cover", "Cover, corrected"),
          "Study 8: bootstrap correction with wrong hyperparameters (logistic outcome, 50 cohorts, 20 refits).",
          "tab:s8b", "tmp4.tex", "llrrrr")
  write_tab(c(L, readLines(file.path(tab_dir, "tmp4.tex"))), "tab_s8_hyper.tex"); file.remove(file.path(tab_dir, "tmp4.tex"))
}

for (f in c("tab_s8", "tab_checks", "tab_s1b", "tab_s7", "tab_study1", "tab_libs", "tab_ist", "tab_s1", "tab_s2", "tab_s3", "tab_s4", "tab_s5", "tab_s6"))
  tryCatch({ get(f)(); cat(f, "ok\n") }, error = function(e) message(f, ": ", conditionMessage(e)))
