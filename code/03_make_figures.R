# ==============================================================================
# 03_make_figures.R  --  figures for the TreeIPF paper (vector PDF, fonts embedded)
#
# base R graphics, Okabe-Ito colours, text 8-9 pt at the printed size
# (170 mm full width), no titles inside the images.
# usage: Rscript 03_make_figures.R [which]   e.g. 1,2,3
# ==============================================================================

suppressPackageStartupMessages(library(data.table))
root_dir <- normalizePath(Sys.getenv("TREEIPF_ROOT", "."), winslash = "/")
source(file.path(root_dir, "code", "00_tree_ipf_lib.R"))
res <- file.path(root_dir, "results"); fig <- file.path(root_dir, "figures")
dir.create(fig, showWarnings = FALSE)
args <- commandArgs(trailingOnly = TRUE)
which_figs <- if (length(args)) as.integer(strsplit(args[1], ",")[[1]]) else 1:6

OI <- c(black = "#000000", orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#999933", blue = "#0072B2", red = "#D55E00", purple = "#CC79A7", grey = "#8C8C8C")
W <- 6.69   # 170 mm

open_pdf <- function(name, h, w = W) {
  cairo_pdf(file.path(fig, name), width = w, height = h, family = "Arial", pointsize = 8.5)
  par(mgp = c(1.7, 0.45, 0), tcl = -0.25, las = 1, cex = 1, lwd = 0.8, bty = "n")
}
panel_tag <- function(t) mtext(t, side = 3, adj = 0, line = 0.4, font = 2, cex = 1.1, at = grconvertX(0.01, "nfc", "user"))
light <- function(col, a = 0.25) adjustcolor(col, alpha.f = a)
axis_clean <- function(side, at = NULL, labels = TRUE, ...) axis(side, at = at, labels = labels, lwd = 0.6, lwd.ticks = 0.6, cex.axis = 1, ...)

# ------------------------------------------------------------------------------
# Figure 1: toy example
# ------------------------------------------------------------------------------
fig1 <- function() {
  set.seed(7)
  n <- 2000; X <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, .6, .6, 1), 2)); colnames(X) <- c("X1", "X2")
  y <- 1 + X[, 1] + 0.6 * X[, 2] + rnorm(n)
  set.seed(9); gb <- fit_cv_gbt(X, y, "gaussian", B_max = 120L, lr = 0.1, maxdepth = 2L, minbucket = 20L, nfolds = 3L)
  set.seed(10); fg <- tree_ipf(gb, M = 20000L, alpha = 0.05)
  lim <- c(-3, 3)
  both <- which(vapply(gb$trees, function(t) all(1:2 %in% t$feat), TRUE))[1:2]
  boxes <- function(tr, r = 1L, b = c(-3, 3, -3, 3)) {
    if (tr$feat[r] == 0L) return(data.frame(x0 = b[1], x1 = b[2], y0 = b[3], y1 = b[4], share = tr$size[r] / tr$size[1]))
    j <- tr$feat[r]; t <- tr$thr[r]
    bl <- b; br <- b
    if (j == 1L) { bl[2] <- min(b[2], t); br[1] <- max(b[1], t) } else { bl[4] <- min(b[4], t); br[3] <- max(b[3], t) }
    rbind(boxes(tr, tr$left[r], bl), boxes(tr, tr$right[r], br))
  }
  open_pdf("fig1_toy_example.pdf", 2.05)
  layout(matrix(1:4, 1), widths = c(1, 1, 1, 1.25)); par(cex = 1)
  for (k in 1:2) {
    tr <- gb$trees[[both[k]]]; col <- c(OI["blue"], OI["red"])[k]
    par(mar = c(2.8, 2.8, 1.6, 0.6))
    plot(X, pch = 16, cex = 0.25, col = light(OI["grey"], 0.35), xlim = lim, ylim = lim, axes = FALSE,
         xlab = expression(X[1]), ylab = expression(X[2]))
    axis_clean(1, at = c(-2, 0, 2)); axis_clean(2, at = c(-2, 0, 2))
    bx <- boxes(tr)
    rect(bx$x0, bx$y0, bx$x1, bx$y1, border = col, lwd = 1.3)
    big <- (pmin(bx$x1, 3) - pmax(bx$x0, -3)) > 0.9 & (pmin(bx$y1, 3) - pmax(bx$y0, -3)) > 0.6
    text(((pmax(bx$x0, -3) + pmin(bx$x1, 3)) / 2)[big], ((pmax(bx$y0, -3) + pmin(bx$y1, 3)) / 2)[big],
         sprintf("%.0f%%", 100 * bx$share)[big], col = col, font = 2)
    mtext(sprintf("tree %d", both[k]), side = 3, line = 0.2, col = col, cex = 1)
    if (k == 1) panel_tag("a")
  }
  # b: points drawn from the distribution rebuilt from all trees
  par(mar = c(2.8, 2.8, 1.6, 0.6))
  idx <- sample.int(nrow(fg$Z), n, TRUE, prob = fg$q)
  plot(fg$Z[idx, ] + matrix(runif(2 * n, -0.04, 0.04), n), pch = 16, cex = 0.25, col = light(OI["blue"], 0.45), xlim = lim, ylim = lim, axes = FALSE,
       xlab = expression(X[1]), ylab = expression(X[2]))
  axis_clean(1, at = c(-2, 0, 2)); axis_clean(2, at = c(-2, 0, 2))
  ins <- fg$cuts
  rect(ins$lo[1], ins$lo[2], ins$hi[1], ins$hi[2], lty = 2, lwd = 1.1, border = OI["black"])
  mtext(sprintf("r: rebuilt %.2f, data %.2f", fg$moments$cor[1, 2], cor(X)[1, 2]), side = 3, line = 0.2, cex = 1)
  panel_tag("b")
  # c: projection inside the outer splits
  par(mar = c(2.8, 3.0, 1.6, 0.6))
  xs <- seq(-3, 3, length.out = 241)
  bg <- fg$Z[sample.int(nrow(fg$Z), 400, TRUE, prob = fg$q), ]
  pd <- sapply(xs, function(v) { b <- bg; b[, 1] <- v; mean(forest_link(gb, b)) })
  lo <- fg$cuts$lo[1]; hi <- fg$cuts$hi[1]
  plot(NA, xlim = lim, ylim = range(pd) + c(-0.4, 0.4), axes = FALSE, xlab = expression(X[1]),
       ylab = "prediction (averaged over X2)")
  rect(lim[1] - 1, -10, lo, 10, col = light(OI["grey"], 0.18), border = NA)
  rect(hi, -10, lim[2] + 1, 10, col = light(OI["grey"], 0.18), border = NA)
  axis_clean(1, at = c(-2, 0, 2)); axis_clean(2)
  m0 <- mean(pd[xs >= lo & xs <= hi]); c0 <- mean(xs[xs >= lo & xs <= hi])
  abline(a = m0 - 1 * c0, b = 1, col = OI["grey"], lwd = 1)
  lines(xs, pd, type = "s", col = OI["black"], lwd = 1)
  segments(lo, m0 + fg$est[2] * (lo - c0), hi, m0 + fg$est[2] * (hi - c0), col = OI["blue"], lwd = 2.2)
  text(lim[2] - 0.02, min(pd) - 0.3, "flat beyond\nouter splits", adj = c(1, 0), col = OI["grey"])
  legend("topleft", inset = c(0.01, 0), c("ensemble", "true slope 1.0", sprintf("TreeIPF slope %.2f", fg$est[2])),
         col = c(OI["black"], OI["grey"], OI["blue"]), lwd = c(1, 1, 2.2), bty = "n", seg.len = 1.4)
  panel_tag("c")
  dev.off()
}

# ------------------------------------------------------------------------------
# Figure 2: where the error comes from, and coverage (Study 1 + Study 6)
# ------------------------------------------------------------------------------
fig2 <- function() {
  pf <- fread(file.path(res, "sim_S1_performance.csv"))
  dc <- fread(file.path(res, "sim_S1_decomposition.csv"))
  bs <- if (file.exists(file.path(res, "sim_S6_bootstrap.csv"))) fread(file.path(res, "sim_S6_bootstrap.csv")) else NULL
  cells <- CJ(rho = c(0, 0.4, 0.7), n = c(300L, 600L, 1500L)); cells[, x := .I + (rho == 0.4) * 0.6 + (rho == 0.7) * 1.2]
  terms <- c(X1 = "X1 (correlated)", X2 = "X2 (correlated)", X3 = "X3 (skewed)", X4 = "X4 (binary)")
  fams <- c(gaussian = "linear", binomial = "logistic")
  open_pdf("fig2_error_and_coverage.pdf", 5.9)
  par(oma = c(0, 1.6, 1.6, 0))
  layout(matrix(1:16, 4, byrow = TRUE)); par(cex = 1)
  for (row in 1:4) {
    fam <- names(fams)[(row - 1) %/% 2 + 1]; meas <- if (row %% 2 == 1) "bias" else "cov"
    for (j in seq_along(terms)) {
      tm <- names(terms)[j]
      par(mar = c(if (row == 4) 2.0 else 1.0, if (j == 1) 3.3 else 1.4, if (row %% 2 == 1) 1.5 else 0.6, 0.3))
      if (meas == "bias") {
        d <- merge(cells, dc[family == fam & term == tm], by = c("rho", "n"))
        tg <- pf[family == fam & term == tm & method == "TreeIPF", .(rho, n, target)]
        d <- merge(d, tg, by = c("rho", "n"))
        ens <- 100 * d$ensemble_part / d$target; tot <- 100 * (d$ensemble_part + d$reconstruction_part) / d$target
        yl <- c(-60, 10)
        plot(NA, xlim = c(0.5, 10.7), ylim = yl, axes = FALSE, xlab = "", ylab = if (j == 1) "relative bias (%)" else "")
        abline(h = 0, col = OI["grey"], lwd = 0.7)
        segments(d$x, 0, d$x, ens, col = light(OI["orange"], 0.9), lwd = 4, lend = 1)
        segments(d$x, ens, d$x, tot, col = OI["blue"], lwd = 4, lend = 1)
        axis_clean(2, at = c(-60, -40, -20, 0), labels = j == 1)
        if (row == 1) mtext(terms[j], side = 3, line = 0.4, cex = 1)
        if (j == 1) mtext(fams[fam], side = 2, line = 4.1, cex = 1, font = 2, las = 0, adj = -0.6)
      } else {
        d <- merge(cells, pf[family == fam & term == tm], by = c("rho", "n"))
        plot(NA, xlim = c(0.5, 10.7), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = if (j == 1) "coverage of 95% CI" else "")
        rect(0, 0.95 - 2 * 0.0154, 11, 0.95 + 2 * 0.0154, col = light(OI["grey"], 0.25), border = NA)
        abline(h = 0.95, col = OI["grey"], lwd = 0.7)
        for (m in list(list("Oracle-interior", OI["black"], 1), list("TreeIPF", OI["blue"], 16))) {
          e <- d[method == m[[1]]][order(x)]
          points(e$x, e$coverage, pch = m[[3]], col = m[[2]], cex = 0.9)
        }
        if (!is.null(bs)) {
          b <- merge(cells, bs[family == fam & term == tm & scenario == "linear"], by = c("rho", "n"))
          points(b$x, b$cov_bc, pch = 17, col = OI["red"], cex = 0.95)
        }
        axis_clean(2, at = c(0, 0.5, 0.95), labels = if (j == 1) c("0", "0.5", "0.95") else FALSE)
      }
      if (row == 4) {
        axis_clean(1, at = cells$x, labels = FALSE)
        mtext(c("\u03c1 = 0", "\u03c1 = 0.4", "\u03c1 = 0.7"), side = 1, line = 0.9, at = tapply(cells$x, cells$rho, mean), cex = 0.95)
      }
    }
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), new = TRUE, mar = c(0, 0, 0, 0)); plot.new()
  legend("top", inset = c(0, 0), horiz = TRUE, bty = "n", x.intersp = 0.6, text.width = NA, xpd = NA,
         legend = c("ensemble shrinkage", "added by rebuilding", "TreeIPF", "TreeIPF, bootstrap-corrected", "oracle inside region"),
         pch = c(15, 15, 16, 17, 1), col = c(OI["orange"], OI["blue"], OI["blue"], OI["red"], OI["black"]), pt.cex = c(1.6, 1.6, 1, 1, 1))
  dev.off()
}

# ------------------------------------------------------------------------------
# Figure 3: rivals -- Type I error, power, accuracy
# ------------------------------------------------------------------------------
METH <- data.table(
  method = c("A:uniform|full|grid", "A:uniform|trim|rootn", "A:marginal|trim|rootn", "PD", "ALE", "TreeIPF", "TreeIPF-Bag", "Oracle"),
  label = c("Naive-Unif", "Unif-rootn", "Marg-IPF", "PD", "ALE", "TreeIPF", "TreeIPF-Bag", "Oracle"),
  col = c(OI["orange"], OI["yellow"], OI["sky"], OI["purple"], OI["green"], OI["blue"], OI["grey"], OI["black"]),
  pch = c(4, 3, 6, 5, 2, 16, 1, 15))

sqrt_axis <- function(at) axis_clean(2, at = sqrt(at), labels = format(at, drop0trailing = TRUE))

fig3 <- function() {
  pf <- fread(file.path(res, "sim_S1_performance.csv"))
  p4 <- fread(file.path(res, "sim_S4_performance.csv"))
  open_pdf("fig3_rivals.pdf", 4.9)
  par(oma = c(1.4, 0, 0, 0))
  layout(matrix(c(1, 1, 1, 2, 3, 4), 2, byrow = TRUE), heights = c(1, 1.05)); par(cex = 1)
  # a: type I error on X5, 18 settings each
  par(mar = c(2.2, 3.4, 1.6, 0.5))
  d <- merge(pf[term == "X5"], METH, by = "method")
  plot(NA, xlim = c(0.5, nrow(METH) + 0.5), ylim = sqrt(c(0, 0.9)), axes = FALSE, xlab = "",
       ylab = "Type I error (X5)")
  rect(0, sqrt(0.02), 99, sqrt(0.08), col = light(OI["grey"], 0.25), border = NA)
  abline(h = sqrt(0.05), col = OI["grey"], lwd = 0.7)
  sqrt_axis(c(0, 0.01, 0.05, 0.2, 0.5, 0.8))
  for (k in seq_len(nrow(METH))) {
    e <- d[method == METH$method[k]]
    xx <- k + ifelse(e$family == "gaussian", -0.15, 0.15) + runif(nrow(e), -0.07, 0.07)
    points(xx, sqrt(e$reject), pch = ifelse(e$family == "gaussian", 16, 1), col = METH$col[k], cex = 0.8)
  }
  axis_clean(1, at = seq_len(nrow(METH)), labels = METH$label, padj = -0.5)
  legend("topright", c("linear", "logistic"), pch = c(16, 1), bty = "n", horiz = TRUE)
  panel_tag("a")
  # b, c: power curves
  pm <- METH[label %in% c("Oracle", "TreeIPF", "Unif-rootn", "Marg-IPF", "Naive-Unif", "PD")]
  for (fam in c("gaussian", "binomial")) {
    par(mar = c(2.8, 3.4, 1.6, 0.5))
    e <- merge(p4[term == "X5" & family == fam], pm, by = "method")
    plot(NA, xlim = range(e$b5), ylim = c(0, 1), axes = FALSE, xlab = expression(beta[5]),
         ylab = if (fam == "gaussian") "rejection rate" else "")
    abline(h = 0.05, col = OI["grey"], lwd = 0.7)
    axis_clean(1, at = sort(unique(e$b5))); axis_clean(2, at = c(0, 0.25, 0.5, 0.75, 1))
    for (k in seq_len(nrow(pm))) {
      f <- e[method == pm$method[k]][order(b5)]
      lines(f$b5, f$reject, col = pm$col[k], lwd = if (pm$label[k] == "TreeIPF") 1.8 else 1)
      points(f$b5, f$reject, col = pm$col[k], pch = pm$pch[k], cex = 0.75)
    }
    mtext(if (fam == "gaussian") "linear" else "logistic", side = 3, line = 0.2, cex = 1)
    if (fam == "gaussian") panel_tag("b")
  }
  # d: RMSE relative to oracle, active coefficients, n = 1500
  par(mar = c(2.8, 6.0, 1.6, 0.5))
  a <- pf[term %in% c("X1", "X2", "X3", "X4") & n == 1500]
  o <- a[method == "Oracle", .(family, rho, term, r0 = rmse)]
  a <- merge(merge(a, o, by = c("family", "rho", "term")), METH, by = "method")
  a <- a[label != "Oracle"]
  lab <- METH[label != "Oracle"]$label
  plot(NA, ylim = c(0.5, length(lab) + 0.5), xlim = log(c(0.8, 12)), axes = FALSE, ylab = "",
       xlab = "RMSE / oracle RMSE")
  abline(v = 0, col = OI["grey"], lwd = 0.7)
  axis_clean(1, at = log(c(1, 2, 4, 8)), labels = c(1, 2, 4, 8))
  axis_clean(2, at = seq_along(lab), labels = lab, tick = FALSE)
  for (k in seq_along(lab)) {
    e <- a[label == lab[k]]
    points(log(e$rmse / e$r0), k + runif(nrow(e), -0.18, 0.18), col = light(e$col[1], 0.8),
           pch = ifelse(e$family == "gaussian", 16, 1), cex = 0.7)
    points(log(median(e$rmse / e$r0)), k, pch = "|", cex = 1.6)
  }
  panel_tag("c")
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), new = TRUE, mar = c(0, 0, 0, 0)); plot.new()
  legend("bottom", legend = pm$label, col = pm$col, pch = pm$pch, lwd = 1, horiz = TRUE, bty = "n",
         x.intersp = 0.5, seg.len = 1.2, xpd = NA)
  dev.off()
}

# ------------------------------------------------------------------------------
# Figure 4: what the files disclose
# ------------------------------------------------------------------------------
fig4 <- function() {
  pv <- fread(file.path(res, "privacy_leakage.csv"))
  s <- pv[, .(L_mean = mean(L_mean), cor_r = mean(cor_r)),
          by = .(cohort, release, depth, minleaf)]
  dcol <- c(`1` = "#9ECAE1", `2` = "#6BAED6", `3` = "#3182BD", `4` = "#08519C", `6` = "#08306B")
  dlty <- c(`1` = 3, `2` = 2, `3` = 1, `4` = 1, `6` = 1)
  open_pdf("fig4_disclosure.pdf", 3.0)
  par(oma = c(1.5, 0, 0, 0))
  layout(matrix(1:4, 1)); par(cex = 1)
  k <- 0
  for (meas in c("cor_r", "L_mean")) for (co in c("Pima", "PhysioNet")) {
    k <- k + 1
    par(mar = c(2.8, if (k %in% c(1, 3)) 3.3 else 1.6, 1.6, 0.4))
    e <- s[cohort == co]
    plot(NA, xlim = log(c(5, 200)), ylim = c(-0.3, 1), axes = FALSE, xlab = "minimum leaf size",
         ylab = if (k == 1) "correlation r" else if (k == 3) "disclosure score" else "")
    abline(h = 0, col = OI["grey"], lwd = 0.7)
    axis_clean(1, at = log(c(5, 10, 25, 50, 100, 200)), labels = c(5, 10, 25, 50, 100, 200))
    axis_clean(2, at = c(0, 0.5, 1), labels = k %in% c(1, 3))
    for (d in names(dcol)) {
      f <- e[release == "leaf counts kept" & depth == as.integer(d)][order(minleaf)]
      lines(log(f$minleaf), f[[meas]], col = dcol[d], lwd = 1.5, lty = dlty[d]); points(log(f$minleaf), f[[meas]], col = dcol[d], pch = 16, cex = 0.6)
    }
    f <- e[release == "counts stripped", .(v = mean(get(meas))), by = minleaf][order(minleaf)]
    lines(log(f$minleaf), f$v, col = OI["red"], lty = 2, lwd = 1.2)
    mtext(sprintf("%s, %s", co, if (meas == "cor_r") "correlations" else "means"), side = 3, line = 0.2, cex = 0.95)
    if (k == 1) panel_tag("a"); if (k == 3) panel_tag("b")
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), new = TRUE, mar = c(0, 0, 0, 0)); plot.new()
  legend("bottom", horiz = TRUE, bty = "n", x.intersp = 0.4, seg.len = 1.6, xpd = NA, text.width = 0.105,
         legend = c(paste("depth", names(dcol)), "counts removed"),
         col = c(dcol, OI["red"]), lwd = c(rep(1.5, 5), 1.2), lty = c(dlty, 2))
  dev.off()
}

# ------------------------------------------------------------------------------
# Figure 5: IST forest plot
# ------------------------------------------------------------------------------
fig5 <- function() {
  ma <- fread(file.path(res, "ist_meta_analysis.csv"))
  lab <- c(age10 = "Age (per 10 years)", female = "Female sex", sbp10 = "Systolic BP (per 10 mmHg)",
           consc_drowsy = "Impaired consciousness", ratrial = "Atrial fibrillation", ct_inf = "Infarct visible on CT",
           rdef1 = "Facial deficit", rdef2 = "Arm or hand deficit")
  st <- data.table(strategy = c("Gold: GLM in all 12", "All 12 from ensembles: TreeIPF",
                                "All 12 from ensembles: TreeIPF + bootstrap", "All 12 from ensembles: Unif-rootn",
                                "All 12 from ensembles: Naive-Unif"),
                   short = c("full-data regressions", "TreeIPF", "TreeIPF + bootstrap", "Unif-rootn", "Naive-Unif"),
                   col = c(OI["black"], OI["blue"], OI["red"], OI["yellow"], OI["orange"]), pch = c(15, 16, 17, 3, 4))
  ord <- names(lab)[order(-ma[strategy == st$strategy[1]][match(names(lab), term)]$OR)]
  open_pdf("fig5_ist_forest.pdf", 4.3)
  layout(matrix(1:2, 1), widths = c(3.3, 1)); par(cex = 1)
  par(mar = c(2.8, 10.5, 0.8, 0.3))
  nS <- nrow(st); gap <- 0.14
  yy <- rev(seq_along(ord))
  plot(NA, xlim = log(c(0.7, 7)), ylim = c(0.5, length(ord) + 0.5), axes = FALSE, xlab = "pooled odds ratio (log scale)", ylab = "")
  abline(v = 0, col = OI["grey"], lwd = 0.7)
  for (i in seq_along(ord)) if (i %% 2 == 0) rect(log(0.5), yy[i] - 0.5, log(10), yy[i] + 0.5, col = light(OI["grey"], 0.1), border = NA)
  axis_clean(1, at = log(c(0.75, 1, 1.5, 2, 3, 5)), labels = c(0.75, 1, 1.5, 2, 3, 5))
  axis_clean(2, at = yy, labels = lab[ord], tick = FALSE)
  for (i in seq_along(ord)) for (k in seq_len(nS)) {
    e <- ma[term == ord[i] & strategy == st$strategy[k]]
    y0 <- yy[i] + (nS + 1) / 2 * gap - k * gap
    segments(log(e$OR_lb), y0, log(e$OR_ub), y0, col = st$col[k], lwd = 1.1)
    points(log(e$OR), y0, pch = st$pch[k], col = st$col[k], cex = 0.8)
  }
  legend("bottomright", st$short, col = st$col, pch = st$pch, lwd = 1, bty = "n", seg.len = 1.2, cex = 0.95, bg = "white")
  # I^2 column
  par(mar = c(2.8, 0.2, 0.8, 0.4))
  plot(NA, xlim = c(0, 100), ylim = c(0.5, length(ord) + 0.5), axes = FALSE, xlab = expression(I^2 ~ "(%)"), ylab = "")
  axis_clean(1, at = c(0, 50, 100))
  for (i in seq_along(ord)) for (k in seq_len(nS)) {
    e <- ma[term == ord[i] & strategy == st$strategy[k]]
    y0 <- yy[i] + (nS + 1) / 2 * gap - k * gap
    segments(0, y0, e$I2, y0, col = light(st$col[k], 0.9), lwd = 2.2, lend = 1)
  }
  dev.off()
}

# ------------------------------------------------------------------------------
# Figure 6: clinical audits over 50 seeds (standardised coefficients)
# ------------------------------------------------------------------------------
fig6 <- function() {
  cf <- fread(file.path(res, "clinical_audit_coefficients.csv"))
  mo <- fread(file.path(res, "clinical_audit_moments.csv"))
  sdv <- mo[seed == 1, .(cohort, term, emp_sd, var_type)]
  cf <- merge(cf[term != "(Intercept)"], sdv, by = c("cohort", "term"))
  cf[, std := ifelse(var_type == "binary", est, est * emp_sd)]
  ms <- data.table(method = c("glm", "TreeIPF", "Unif-rootn", "Naive-Unif"),
                   col = c(OI["black"], OI["blue"], OI["yellow"], OI["orange"]), pch = c(18, 16, 3, 4))
  nice <- c(glu = "Glucose", bmi = "Body mass index", ped = "Pedigree function", npreg = "Pregnancies", age = "Age",
            skin = "Skinfold", bp = "Diastolic BP", age10 = "Age", female = "Female sex", hr10 = "Highest heart rate",
            sbp10 = "Lowest systolic BP", temp_max = "Highest temperature", gcs_min = "Lowest GCS", bun10 = "Highest urea",
            creat_max = "Highest creatinine", wbc10 = "Highest white cell count", mechvent = "Ventilation")
  open_pdf("fig6_clinical_audits.pdf", 3.5)
  layout(matrix(1:2, 1), widths = c(1, 1.15)); par(cex = 1)
  for (co in c("Pima", "PhysioNet")) {
    e <- cf[cohort == co]
    g <- e[method == "glm", .(std = mean(std)), by = term][order(std)]
    tl <- g$term
    par(mar = c(2.8, 8.6, 1.6, 0.6))
    plot(NA, xlim = range(e[method %in% ms$method]$std), ylim = c(0.5, length(tl) + 0.5), axes = FALSE,
         xlab = "standardised log odds ratio", ylab = "")
    abline(v = 0, col = OI["grey"], lwd = 0.7)
    axis_clean(1); axis_clean(2, at = seq_along(tl), labels = nice[tl], tick = FALSE)
    for (k in seq_len(nrow(ms))) {
      f <- e[method == ms$method[k], .(lo = quantile(std, 0.05), mid = median(std), hi = quantile(std, 0.95)), by = term]
      y0 <- match(f$term, tl) + (2.5 - k) * 0.17
      segments(f$lo, y0, f$hi, y0, col = ms$col[k], lwd = 1.1)
      points(f$mid, y0, pch = ms$pch[k], col = ms$col[k], cex = if (k == 1) 1.2 else 0.8)
    }
    mtext(co, side = 3, line = 0.2, cex = 1)
    panel_tag(if (co == "Pima") "a" else "b")
  }
  legend("bottomright", c("full-data regression", "TreeIPF", "Unif-rootn", "Naive-Unif"), col = ms$col, pch = ms$pch,
         lwd = 1, bty = "n", seg.len = 1)
  dev.off()
}

for (k in which_figs) {
  f <- get(paste0("fig", k))
  ok <- tryCatch({ f(); TRUE }, error = function(e) { message("figure ", k, ": ", conditionMessage(e)); FALSE })
  if (ok) cat("figure", k, "written\n")
}
