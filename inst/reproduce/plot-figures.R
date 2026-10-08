# Draws the published figures from the values that inst/reproduce/reproduce-figures.R
# writes to inst/extdata/published-figures/. Each function reads its CSV files from the
# folder dir and draws one figure on the current graphics device, with the axis ranges,
# tick marks, line types and symbols of the original figure as far as base graphics allow.
# The vignette 'Validation by reproducing published figures' sources this file;
# write_figure_pngs() writes all figures as PNG files.

fig_dir <- function() system.file('extdata', 'published-figures', package = 'bbssr')
read_fig <- function(dir, name) utils::read.csv(file.path(dir, name))
# Axis with labelled ticks at 'at' and unlabelled minor ticks at 'minor'
fig_axis <- function(side, at, labels = TRUE, minor = NULL, ...) {
  graphics::axis(side, at = at, labels = labels, ...)
  if (!is.null(minor)) graphics::axis(side, at = minor, labels = FALSE, tcl = -0.25)
}
# Grey strip above a panel, as in the facets of the figures of Kieser (2020)
fig_strip <- function(label) {
  usr <- graphics::par('usr')
  h <- 0.09 * (usr[4] - usr[3])
  graphics::rect(usr[1], usr[4], usr[2], usr[4] + h, col = 'grey85', border = 'grey40',
                 xpd = NA)
  graphics::text(mean(usr[1:2]), usr[4] + h / 2, label, xpd = NA, cex = 0.9)
}
# Label such as 'a)' at the left edge of the device, level with the top of the panel
fig_tag <- function(tag, above = 0.04) {
  usr <- graphics::par('usr')
  graphics::text(graphics::grconvertX(0.01, 'ndc', 'user'),
                 usr[4] + above * (usr[4] - usr[3]), tag, adj = c(0, 0.5), cex = 1.1,
                 xpd = NA)
}
# Axis limits of ggplot2: the range of the data widened by 5 per cent on each side
fig_expand <- function(x) range(x) + c(-1, 1) * 0.05 * diff(range(x))
fig_restore <- function() {
  op <- graphics::par(no.readonly = TRUE)
  function() {
    graphics::layout(1)
    graphics::par(op)
  }
}
n1_label <- function(n1) as.expression(lapply(n1, function(v) {
  bquote(italic(n)[1 * '.'] == .(v))
}))

# Friede and Kieser (2004) ---------------------------------------------------------------
plot_fk2004_figure1 <- function(dir = fig_dir()) {
  restore <- fig_restore()
  on.exit(restore())
  d <- read_fig(dir, 'friede-kieser-2004-figure1.csv')
  pilot <- seq(20, 100, by = 20)
  col <- grey(c(0.15, 0.35, 0.5, 0.65, 0.8))
  # Panels in cells 1 and 3, their legends below in cells 2 and 4
  graphics::layout(matrix(c(1, 3, 2, 4), 2, byrow = TRUE), heights = c(5, 1.1))
  for (rate in c(0.1, 0.5)) {
    s <- d[d$pi == rate, ]
    graphics::par(mar = c(4, 4.5, 2.5, 1), las = 1)
    plot(NA, xlim = c(0, 200), ylim = c(0, 0.08), xaxs = 'i', yaxs = 'i', axes = FALSE,
         xlab = expression(italic(n)), ylab = expression(alpha[act]),
         main = bquote(pi == .(rate)), font.main = 1)
    fig_axis(1, seq(0, 200, by = 20), minor = seq(0, 200, by = 5))
    fig_axis(2, seq(0, 0.08, by = 0.01), sprintf('%.2f', seq(0, 0.08, by = 0.01)),
             minor = seq(0, 0.08, by = 0.0025))
    graphics::box()
    graphics::abline(h = 0.05, col = 'grey45')
    for (k in rev(seq_along(pilot))) {
      graphics::lines(s$n, s[[paste0('pilot.', pilot[k])]], col = col[k], lwd = 2)
    }
    graphics::lines(s$n, s$fixed, lty = 2, lwd = 2)
    # Legend of two rows under each panel, filled by column
    graphics::par(mar = c(0, 4.5, 0, 1))
    plot.new()
    ord <- c(1, 4, 2, 5, 3, 6)
    graphics::legend('center', ncol = 3, lwd = 2,
                     legend = c('fixed', n1_label(pilot))[ord],
                     lty = c(2, rep(1, 5))[ord], col = c('black', col)[ord])
  }
}

plot_fk2004_figure2 <- function(dir = fig_dir()) {
  restore <- fig_restore()
  on.exit(restore())
  d <- read_fig(dir, 'friede-kieser-2004-figure2.csv')
  graphics::par(mfrow = c(1, 2), mar = c(4.5, 5.5, 2.5, 1), las = 1)
  for (th in c(1, 3)) {
    s <- d[d$theta == th, ]
    plot(NA, xlim = c(20, 200), ylim = c(0, 0.07), xaxs = 'i', yaxs = 'i', axes = FALSE,
         xlab = expression(italic(n)[1 * '.']), ylab = '', main = bquote(theta == .(th)),
         font.main = 1)
    graphics::title(ylab = expression(alpha[act]), line = 4)
    fig_axis(1, seq(20, 200, by = 20), minor = seq(20, 200, by = 10))
    fig_axis(2, seq(0, 0.07, by = 0.01), sprintf('%.4f', seq(0, 0.07, by = 0.01)),
             minor = seq(0, 0.07, by = 0.0025))
    graphics::box()
    graphics::abline(h = 0.05, col = 'grey45')
    # As drawn in the article: the fixed design dashed and the internal pilot study design
    # solid (the caption of the article states the reverse)
    for (v in c('fixed.min', 'fixed.mean', 'fixed.max')) {
      graphics::lines(s$pilot, s[[v]], type = 'o', lty = 2, pch = 20, lwd = 2)
    }
    for (v in c('ips.min', 'ips.mean', 'ips.max')) {
      graphics::lines(s$pilot, s[[v]], type = 'o', lty = 1, pch = 20, lwd = 2)
    }
  }
}

plot_fk2004_figure3 <- function(dir = fig_dir()) {
  restore <- fig_restore()
  on.exit(restore())
  d <- read_fig(dir, 'friede-kieser-2004-figure3.csv')
  pilot <- c(40, 80, 120)
  col <- grey(c(0, 0.5, 0.75))
  graphics::layout(matrix(c(1, 3, 2, 4), 2, byrow = TRUE), heights = c(5, 1.1))
  for (th in c(1, 3)) {
    s <- d[d$theta == th, ]
    xmax <- if (th == 1) 0.5 else 0.9
    graphics::par(mar = c(4, 4.5, 2.5, 1), las = 1)
    plot(NA, xlim = c(0, xmax), ylim = c(0.68, 0.90), xaxs = 'i', yaxs = 'i',
         axes = FALSE, xlab = expression(pi[1]), ylab = 'Power',
         main = bquote(theta == .(th)), font.main = 1)
    fig_axis(1, seq(0, xmax, by = 0.1), sprintf('%.1f', seq(0, xmax, by = 0.1)),
             minor = seq(0, xmax, by = 0.025))
    fig_axis(2, seq(0.68, 0.90, by = 0.02), sprintf('%.2f', seq(0.68, 0.90, by = 0.02)),
             minor = seq(0.68, 0.90, by = 0.005))
    graphics::box()
    graphics::abline(h = 0.8, col = 'grey45')
    for (k in rev(seq_along(pilot))) {
      graphics::lines(s$pi1, s[[paste0('pilot.', pilot[k])]], type = 'o', pch = 20,
                      lwd = 2.5, col = col[k])
    }
    graphics::lines(s$pi1, s$fixed, type = 'o', lty = 2, pch = 20, lwd = 2.5)
    graphics::par(mar = c(0, 4.5, 0, 1))
    plot.new()
    ord <- c(1, 3, 2, 4)
    graphics::legend('center', ncol = 2, lwd = 2.5, pch = 20,
                     legend = c('fixed', n1_label(pilot))[ord],
                     lty = c(2, 1, 1, 1)[ord], col = c('black', col)[ord])
  }
}

plot_fk2004_figure4 <- function(dir = fig_dir()) {
  restore <- fig_restore()
  on.exit(restore())
  a <- read_fig(dir, 'friede-kieser-2004-figure4-fixed.csv')
  b <- read_fig(dir, 'friede-kieser-2004-figure4-recalculation.csv')
  graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.8, 1, 1), las = 1)
  plot(NA, xlim = c(0, 0.5), ylim = c(0, 0.06), xaxs = 'i', yaxs = 'i', axes = FALSE,
       xlab = expression(pi), ylab = expression(alpha[act]^fix))
  fig_axis(1, seq(0, 0.5, by = 0.05), sprintf('%.2f', seq(0, 0.5, by = 0.05)),
           minor = seq(0, 0.5, by = 0.0125), cex.axis = 0.75)
  fig_axis(2, seq(0, 0.06, by = 0.01), sprintf('%.2f', seq(0, 0.06, by = 0.01)),
           minor = seq(0, 0.06, by = 0.0025))
  graphics::box()
  graphics::abline(h = 0.05, col = 'grey45')
  graphics::lines(a$pi, a$level, lwd = 2.5)
  plot(NA, xlim = c(0, 400), ylim = c(0, 0.06), xaxs = 'i', yaxs = 'i', axes = FALSE,
       xlab = expression(italic(n)), ylab = expression(alpha[act]^recalc))
  fig_axis(1, seq(0, 400, by = 100), minor = seq(0, 400, by = 25))
  fig_axis(2, seq(0, 0.06, by = 0.01), sprintf('%.2f', seq(0, 0.06, by = 0.01)),
           minor = seq(0, 0.06, by = 0.0025))
  graphics::box()
  graphics::abline(h = 0.05, col = 'grey45')
  graphics::lines(b$n, b$level, lwd = 2.5)
}

# Kieser (2020) ------------------------------------------------------------------------
# Panels of Figure 21.1. The axis limits are read from the book, and the dotted lines are
# drawn where the book draws them, at 0.02475 and 0.0275 (its caption gives 0.0225 and
# 0.0275). 'stop' adds the levels under the rule that stops the trial after the pilot
kieser_scatter <- function(dir, stop) {
  restore <- fig_restore()
  on.exit(restore())
  d <- read_fig(dir, 'kieser-2020-figure21-1.csv')
  lims <- list('1' = c(0.0194, 0.0316), '3' = c(0.0163, 0.0317))
  graphics::par(mfrow = c(2, 3), mar = c(3.4, 3.4, 2.2, 0.5), oma = c(0, 1.6, 0, 0),
                las = 1, mgp = c(2, 0.6, 0))
  for (ri in c(1, 3)) {
    for (fr in c(0.25, 0.5, 0.75)) {
      s <- d[d$r == ri & d$fraction == fr, ]
      lim <- lims[[as.character(ri)]]
      plot(NA, xlim = lim, ylim = lim, xaxs = 'i', yaxs = 'i', axes = FALSE, xlab = '',
           ylab = '')
      for (side in 1:2) {
        graphics::axis(side, at = c(0.020, 0.024, 0.028),
                       labels = c('0.020', '0.024', '0.028'), cex.axis = 0.8,
                       col = 'grey40')
      }
      graphics::box(col = 'grey40')
      graphics::abline(0, 1, lty = 2)
      graphics::abline(v = c(0.02475, 0.0275), h = c(0.02475, 0.0275), lty = 3,
                       col = 'grey30')
      if (stop) {
        graphics::points(s$fixed, s$ips, pch = 19, cex = 0.8, col = 'grey70')
        graphics::points(s$fixed, s$ips.stop, pch = 19, cex = 0.8)
      } else {
        graphics::points(s$fixed, s$ips, pch = 19, cex = 0.8)
      }
      fig_strip(format(fr))
      # Axis titles once per row, as in the book
      if (fr == 0.25) {
        fig_tag(if (ri == 1) 'a)' else 'b)', above = 0.045)
        graphics::mtext('actual level recalculation design', side = 2, line = 2.6,
                        las = 0, cex = 0.7)
      }
      if (fr == 0.5) {
        graphics::mtext('actual level fixed design', side = 1, line = 2.1, cex = 0.75)
      }
    }
  }
}
plot_kieser_figure21_1 <- function(dir = fig_dir()) kieser_scatter(dir, stop = FALSE)
plot_kieser_figure21_1_stop <- function(dir = fig_dir()) kieser_scatter(dir, stop = TRUE)

# Two stacked panels with the legend 'design' on the right, as in Figures 21.2 and 21.3
kieser_curves <- function(panels, ylab, hline) {
  restore <- fig_restore()
  on.exit(restore())
  graphics::layout(matrix(c(1, 2, 3, 3), 2), widths = c(5, 1))
  for (pn in panels) {
    graphics::par(mar = c(4, 5, 1.5, 0.5), las = 1)
    plot(NA, xlim = fig_expand(pn$x), ylim = fig_expand(c(pn$y1, pn$y2, hline)),
         xaxs = 'i', yaxs = 'i', axes = FALSE, xlab = 'p', ylab = '')
    graphics::title(ylab = ylab, line = 3.8)
    graphics::axis(1, at = pn$xat, labels = pn$xlab, col = 'grey40', cex.axis = 0.85)
    graphics::axis(2, at = pn$yat, labels = pn$ylab, col = 'grey40', cex.axis = 0.85)
    graphics::box(col = 'grey40')
    graphics::abline(h = hline, col = 'grey60')
    graphics::lines(pn$x, pn$y1, col = 'grey20', lwd = 1.5)
    graphics::lines(pn$x, pn$y2, col = 'grey20', lwd = 1.5, lty = 3)
    fig_tag(pn$tag)
  }
  graphics::par(mar = c(0, 0, 0, 0))
  plot.new()
  graphics::legend('center', title = 'design', legend = c('fixed', 'recalculation'),
                   lty = c(1, 3), lwd = 1.5, col = 'grey20', bty = 'n', cex = 0.9)
}

plot_kieser_figure21_2 <- function(dir = fig_dir()) {
  d <- read_fig(dir, 'kieser-2020-figure21-2.csv')
  d <- d[d$shown, ]
  a <- d[d$Delta == -0.15, ]
  b <- d[d$Delta == -0.30, ]
  kieser_curves(list(
    list(x = a$p, y1 = a$fixed, y2 = a$recalculation, tag = 'a)',
         xat = c(0.25, 0.5, 0.75), xlab = c('0.25', '0.50', '0.75'),
         yat = c(0.024, 0.025, 0.026), ylab = c('0.024', '0.025', '0.026')),
    list(x = b$p, y1 = b$fixed, y2 = b$recalculation, tag = 'b)',
         xat = c(0.2, 0.4, 0.6, 0.8), xlab = c('0.2', '0.4', '0.6', '0.8'),
         yat = seq(0.025, 0.027, by = 0.0005),
         ylab = sprintf('%.4f', seq(0.025, 0.027, by = 0.0005)))
  ), ylab = 'actual level', hline = 0.025)
}

plot_kieser_figure21_3 <- function(dir = fig_dir()) {
  d <- read_fig(dir, 'kieser-2020-figure21-3.csv')
  panel <- function(ri, tag) {
    s <- d[d$r == ri, ]
    list(x = s$p, y1 = s$fixed, y2 = s$recalculation, tag = tag,
         xat = c(0.25, 0.5, 0.75), xlab = c('0.25', '0.50', '0.75'),
         yat = seq(0.75, 1, by = 0.05), ylab = sprintf('%.2f', seq(0.75, 1, by = 0.05)))
  }
  kieser_curves(list(panel(1, 'a)'), panel(3, 'b)')), ylab = 'power', hline = 0.8)
}

plot_kieser_figure21_4 <- function(dir = fig_dir()) {
  restore <- fig_restore()
  on.exit(restore())
  b <- read_fig(dir, 'kieser-2020-figure21-4.csv')
  o <- read_fig(dir, 'kieser-2020-figure21-4-outside.csv')
  w <- 0.0225
  graphics::par(mar = c(4, 4.5, 0.5, 0.5), las = 1)
  plot(NA, xlim = fig_expand(c(b$p - w, b$p + w)),
       ylim = fig_expand(c(b$whisker.low, b$whisker.high, o$N)), xaxs = 'i', yaxs = 'i',
       axes = FALSE, xlab = 'p', ylab = 'sample size')
  graphics::axis(1, at = c(0.25, 0.5, 0.75), labels = c('0.25', '0.50', '0.75'),
                 col = 'grey40', cex.axis = 0.85)
  graphics::axis(2, at = seq(150, 350, by = 50), col = 'grey40', cex.axis = 0.85)
  graphics::box(col = 'grey40')
  graphics::segments(b$p, b$N.q75, b$p, b$whisker.high)
  graphics::segments(b$p, b$N.q25, b$p, b$whisker.low)
  graphics::rect(b$p - w, b$N.q25, b$p + w, b$N.q75, col = 'white')
  graphics::segments(b$p - w, b$N.q50, b$p + w, b$N.q50, lwd = 3)
  graphics::points(o$p, o$N, pch = 8, cex = 0.45)
  graphics::points(b$p, b$E.N, pch = 4, cex = 1.3)
}

# Friede, Mitchell and Mueller-Velten (2007) ---------------------------------------------
plot_fmm2007_figure1 <- function(dir = fig_dir()) {
  restore <- fig_restore()
  on.exit(restore())
  d <- read_fig(dir, 'friede-2007-figure1.csv')
  graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 1, 1), pty = 's')
  for (tst in c('Blackwelder', 'Farrington-Manning')) {
    s <- d[d$Test == tst, ]
    plot(NA, xlim = c(0.02, 0.03), ylim = c(0.02, 0.03), axes = FALSE,
         xlab = 'Fixed size', ylab = 'Sample size reestimation')
    for (side in 1:2) {
      graphics::axis(side, at = seq(0.02, 0.03, by = 0.002), labels = FALSE)
      graphics::axis(side, at = c(0.020, 0.024, 0.028),
                     labels = c('0.020', '0.024', '0.028'), tick = FALSE)
    }
    graphics::box()
    graphics::abline(v = c(0.0225, 0.0275), h = c(0.0225, 0.0275), lty = 3,
                     col = 'grey40')
    graphics::abline(0, 1, lty = 3, col = 'grey40')
    graphics::points(s$fixed, s$reestimation, pch = 1, cex = 0.7)
  }
}

# Axes of Figures 2 and 3: the label 'Power' above the y axis, as in the article
fmm_frame <- function(ylim, yat, yminor, digits, main = NULL) {
  plot(NA, xlim = c(-0.2, 0.2), ylim = ylim, yaxs = 'i', axes = FALSE, ylab = '',
       xlab = expression(pi - pi^a), main = main, font.main = 1)
  fig_axis(1, seq(-0.2, 0.2, by = 0.1), sprintf('%.1f', seq(-0.2, 0.2, by = 0.1)),
           minor = seq(-0.2, 0.2, by = 0.025))
  fig_axis(2, yat, formatC(yat, format = 'f', digits = digits), minor = yminor)
  graphics::box()
  graphics::mtext('Power', side = 3, at = graphics::par('usr')[1], line = 0.1, cex = 0.8)
}

plot_fmm2007_figure2 <- function(dir = fig_dir()) {
  restore <- fig_restore()
  on.exit(restore())
  d <- read_fig(dir, 'friede-2007-figure2.csv')
  graphics::par(mfrow = c(1, 2), mar = c(4.5, 4, 3.5, 1), las = 1)
  for (q in c(1, 1 / 3)) {
    fmm_frame(c(0.7, 1), seq(0.7, 1, by = 0.05), seq(0.7, 1, by = 0.01), 2,
              main = if (q == 1) expression(theta == 1) else expression(theta == 0.33))
    for (pa in c(0.5, 0.7)) {
      s <- d[abs(d$q - q) < 1e-9 & d$pa == pa, ]
      cl <- if (pa == 0.5) 'black' else 'grey55'
      graphics::lines(s$shift, s$fixed, type = 'o', pch = 19, lty = 2, lwd = 2, col = cl,
                      xpd = TRUE)
      graphics::lines(s$shift, s$reestimation, type = 'o', pch = 19, lwd = 2, col = cl,
                      xpd = TRUE)
    }
  }
}

plot_fmm2007_figure3 <- function(dir = fig_dir()) {
  restore <- fig_restore()
  on.exit(restore())
  d <- read_fig(dir, 'friede-2007-figure3.csv')
  graphics::par(mfrow = c(2, 2), mar = c(4.5, 4, 2.5, 1), oma = c(0, 2, 1.5, 0), las = 1)
  for (pa in c(0.5, 0.7)) {
    for (design in c('fixed', 'reestimation')) {
      fmm_frame(c(0.5, 1), seq(0.5, 1, by = 0.1), seq(0.5, 1, by = 0.02), 1)
      graphics::abline(h = 0.8)
      for (d1 in c(-0.02, -0.01, 0.01, 0.02, 0)) {
        s <- d[d$pa == pa & abs(d$delta1 - d1) < 1e-9, ]
        graphics::lines(s$shift, s[[design]], type = 'o', pch = 19, lwd = 2, xpd = TRUE,
                        lty = if (d1 == 0) 1 else 2,
                        col = if (abs(d1) > 0.015) 'grey55' else 'black')
      }
      if (pa == 0.5) {
        graphics::mtext(if (design == 'fixed') 'Fixed' else 'Reestimation', side = 3,
                        line = 1.8)
      }
      if (design == 'fixed') {
        graphics::mtext(bquote(pi^a == .(pa)), side = 2, line = 3.5, las = 0)
      }
    }
  }
}

# Boschloo (1970) ----------------------------------------------------------------------
plot_boschloo1970_figure2 <- function(dir = fig_dir()) {
  restore <- fig_restore()
  on.exit(restore())
  d <- read_fig(dir, 'boschloo-1970-figure2.csv')
  top <- read_fig(dir, 'boschloo-1970-figure2-maxima.csv')
  graphics::par(mar = c(4.5, 5, 1, 1), las = 1)
  plot(NA, xlim = c(0, 1), ylim = c(0, 0.07), xaxs = 'i', yaxs = 'i', axes = FALSE,
       xlab = 'p', ylab = 'unconditional level of significance')
  graphics::axis(1, at = seq(0, 1, by = 0.1),
                 labels = c(sub('^0', '', sprintf('%.1f', seq(0, 0.9, by = 0.1))), '1'))
  graphics::axis(2, at = seq(0, 0.07, by = 0.01),
                 labels = sub('^0', '', sprintf('%.2f', seq(0, 0.07, by = 0.01))))
  graphics::box()
  graphics::abline(h = 0.05)
  graphics::lines(d$p, d$alpha.gamma, lwd = 2)
  graphics::lines(d$p, d$alpha, lwd = 2)
  for (r in c(6, 10, 17, 21)) graphics::lines(d$p, d[[paste0('f', r)]])
  pos <- top$alpha.r > 0
  graphics::points(top$p[pos], top$f.max[pos], pch = 1, cex = 0.7)
  k <- which.max(d$alpha.gamma)
  graphics::text(d$p[k], d$alpha.gamma[k], expression(alpha[gamma](italic(p))), pos = 1)
  k <- which.max(d$alpha)
  graphics::text(d$p[k], d$alpha[k], expression(alpha(italic(p))), pos = 3)
  for (r in c(6, 10, 17, 21)) {
    graphics::text(r / 25, top$f.max[r + 1], bquote(italic(f)[.(r)](italic(p))), pos = 3,
                   cex = 0.85)
  }
  graphics::legend('topright', bty = 'n', cex = 0.9, legend = expression(
    italic(m) == 15, italic(n) == 10, italic(t) == 25, alpha == '.05', gamma == '.09'))
}

# All figures as PNG files in the folder out.dir
write_figure_pngs <- function(out.dir, dir = fig_dir()) {
  dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
  figs <- list(
    'friede-kieser-2004-figure1' = list(plot_fk2004_figure1, 1800, 1000),
    'friede-kieser-2004-figure2' = list(plot_fk2004_figure2, 1800, 900),
    'friede-kieser-2004-figure3' = list(plot_fk2004_figure3, 1800, 1000),
    'friede-kieser-2004-figure4' = list(plot_fk2004_figure4, 1800, 850),
    'kieser-2020-figure21-1' = list(plot_kieser_figure21_1, 1600, 1150),
    'kieser-2020-figure21-1-pilot-stop' = list(plot_kieser_figure21_1_stop, 1600, 1150),
    'kieser-2020-figure21-2' = list(plot_kieser_figure21_2, 1600, 1350),
    'kieser-2020-figure21-3' = list(plot_kieser_figure21_3, 1600, 1350),
    'kieser-2020-figure21-4' = list(plot_kieser_figure21_4, 1600, 750),
    'friede-2007-figure1' = list(plot_fmm2007_figure1, 1600, 850),
    'friede-2007-figure2' = list(plot_fmm2007_figure2, 1300, 950),
    'friede-2007-figure3' = list(plot_fmm2007_figure3, 1500, 1500),
    'boschloo-1970-figure2' = list(plot_boschloo1970_figure2, 1400, 1000)
  )
  for (nm in names(figs)) {
    grDevices::png(file.path(out.dir, paste0(nm, '.png')), width = figs[[nm]][[2]],
                   height = figs[[nm]][[3]], res = 150)
    figs[[nm]][[1]](dir)
    grDevices::dev.off()
  }
  invisible(file.path(out.dir, paste0(names(figs), '.png')))
}
