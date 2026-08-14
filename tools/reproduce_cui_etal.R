###############################################################################
# reproduce_cui_etal.R
#
# Reproduces the numerical results of
#
#   Cui, Y., Homma, G. and Daimon, T.
#   "Exact unconditional tests for blinded sample size re-estimation"
#
# using bbssr 2.0.0.
#
# Version 2.0.0 corrects the handling of ties in the ordering statistic of the
# exact unconditional tests. Outcomes that share the same value of the ordering
# statistic now receive the same p-value, which changes the rejection region at
# the sample sizes where a tie sits next to the critical value. The figures
# produced here can therefore differ from the published ones. Every number the
# script prints is computed, none is taken from the article, and the article's
# reported values are shown alongside for comparison.
#
# Run from the package root with
#
#   source('tools/reproduce_cui_etal.R')
#
# Run time: the full grid takes a few hours with QUICK set to FALSE. Results are
# cached one setting at a time, so the script can be interrupted and restarted.
###############################################################################

library(bbssr)
library(ggplot2)

## Configuration --------------------------------------------------------------

# TRUE runs a coarse subset that finishes in a few minutes
QUICK <- TRUE

# Directory for cached results, tables and figures
OUT.DIR <- 'reproduce-output'

# Number of grid points for the nuisance parameter of the unconditional tests.
# Version 1 of bbssr used 100, which is what the article was computed with. A
# finer grid finds a larger maximum and therefore gives more conservative tests.
N.GRID <- 100

# Berger-Boos parameter. Set to 1e-4 to reproduce Web Appendix C
BB.GAMMA <- 0

## Settings of the article ----------------------------------------------------
#
# Section 4 fixes the response probability of the experimental group at
# theta.E = 0.4 and varies the treatment effect, so theta.C = 0.4 - Delta. The
# nominal one-sided level is 0.025 and the target power is 0.9. The nuisance
# parameter runs over [r Delta / (1 + r), 1 - Delta / (1 + r)] in steps of 0.01,
# which is the range that keeps both group probabilities inside the unit
# interval when the treatment effect is present. Re-estimation uses the
# unrestricted rule at information fraction omega.

THETA.E <- 0.4
ALPHA <- 0.025
TAR.POWER <- 0.9

if (QUICK) {
  TESTS <- c('Boschloo', 'Z-pool')
  R.VALUES <- c(1, 2)
  DELTAS <- c(0.30, 0.36)
  OMEGAS <- c(0.3, 0.8, 0.9)
  THETA.STEP <- 0.05
} else {
  TESTS <- c('Boschloo', 'Z-pool')
  R.VALUES <- c(1, 2)
  DELTAS <- seq(0.24, 0.36, by = 0.03)
  OMEGAS <- seq(0.1, 0.9, by = 0.1)
  THETA.STEP <- 0.01
}

## Mapping between the article and bbssr --------------------------------------
#
# Group 1 of bbssr is the experimental group E and group 2 is the control group
# C, so r is the allocation ratio to the experimental group in both. The article
# writes the sample size as a total n0 and splits it as nC = n0 / (1 + r) and
# nE = r nC, whereas bbssr carries the two group sizes. With an integer r the
# two agree, and BinarySampleSize() reproduces Algorithm 1 of Web Appendix A
# exactly: the initial value is the same large-sample approximation and the
# search moves the control group by one patient, that is the total by (1 + r).
#
# The stage 1 sizes differ for r > 1. The article sets nC1 = ceiling(omega nC)
# and nE1 = r nC1, while bbssr sets nE1 = ceiling(omega nE) = ceiling(r omega
# nC). These coincide whenever ceiling(r omega nC) equals r ceiling(omega nC),
# which always holds for r = 1 and holds for most but not all omega when r = 2.
# check_stage1_agreement() below reports the settings where they part company,
# and those settings are reproduced only approximately.

theta.grid <- function(Delta, r, step = THETA.STEP) {
  lower <- r * Delta / (1 + r)
  upper <- 1 - Delta / (1 + r)
  round(seq(lower, upper, by = step), 10)
}

## Fixed-sample size ----------------------------------------------------------

message('Computing the fixed-sample size for every setting')

n0.table <- do.call(rbind, lapply(TESTS, function(tst) {
  do.call(rbind, lapply(R.VALUES, function(r) {
    do.call(rbind, lapply(DELTAS, function(Delta) {
      ss <- BinarySampleSize(
        p1 = THETA.E, p2 = THETA.E - Delta, r = r,
        alpha = ALPHA, tar.power = TAR.POWER, Test = tst,
        n.grid = N.GRID, bb.gamma = BB.GAMMA
      )
      data.frame(Test = tst, r = r, Delta = Delta,
                 theta.E = THETA.E, theta.C = THETA.E - Delta,
                 nE0 = ss$N1, nC0 = ss$N2, n0 = ss$N,
                 power.at.n0 = ss$Power, stringsAsFactors = FALSE)
    }))
  }))
}))

## Diagnostic on the stage 1 sizes --------------------------------------------

check_stage1_agreement <- function(n0.table, omegas) {
  out <- do.call(rbind, lapply(seq_len(nrow(n0.table)), function(k) {
    row <- n0.table[k, ]
    data.frame(
      Test = row$Test, r = row$r, Delta = row$Delta, omega = omegas,
      article.nE1 = row$r * ceiling(omegas * row$nC0),
      bbssr.nE1 = ceiling(omegas * row$nE0),
      stringsAsFactors = FALSE
    )
  }))
  out$agree <- out$article.nE1 == out$bbssr.nE1
  out
}

stage1.check <- check_stage1_agreement(n0.table, OMEGAS)

## One setting of the grid ----------------------------------------------------

run_setting <- function(tst, r, Delta, omega, n0.row) {
  theta <- theta.grid(Delta, r)
  common <- list(
    p = theta, Delta.A = Delta,
    N1 = n0.row$nE0, N2 = n0.row$nC0,
    omega = omega, r = r,
    alpha = ALPHA, tar.power = TAR.POWER, Test = tst,
    restricted = FALSE, n.grid = N.GRID, bb.gamma = BB.GAMMA
  )
  # Delta.T = 0 puts both groups at the same response probability, so the two
  # power columns become rejection probabilities under the null hypothesis
  null.run <- do.call(BinaryPowerBSSR, c(common, list(Delta.T = 0)))
  alt.run <- do.call(BinaryPowerBSSR, c(common, list(Delta.T = Delta)))
  data.frame(
    Test = tst, r = r, Delta = Delta, omega = omega,
    theta = null.run$p,
    TIE.BSSR = null.run$power.BSSR,
    TIE.fixed = null.run$power.TRAD,
    power.BSSR = alt.run$power.BSSR,
    power.fixed = alt.run$power.TRAD,
    EN.BSSR = alt.run$E.N,
    n0 = n0.row$n0,
    stringsAsFactors = FALSE
  )
}

## Main loop ------------------------------------------------------------------

dir.create(file.path(OUT.DIR, 'cache'), recursive = TRUE, showWarnings = FALSE)

grid <- expand.grid(omega = OMEGAS, Delta = DELTAS, r = R.VALUES, Test = TESTS,
                    stringsAsFactors = FALSE)

message(sprintf('Running %d settings', nrow(grid)))

results <- vector('list', nrow(grid))
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]
  key <- sprintf('%s_r%g_d%g_w%g_g%d_bb%g',
                 g$Test, g$r, g$Delta, g$omega, N.GRID, BB.GAMMA)
  file <- file.path(OUT.DIR, 'cache', paste0(key, '.rds'))
  if (file.exists(file)) {
    results[[k]] <- readRDS(file)
    next
  }
  n0.row <- n0.table[n0.table$Test == g$Test & n0.table$r == g$r &
                       n0.table$Delta == g$Delta, ]
  started <- Sys.time()
  results[[k]] <- run_setting(g$Test, g$r, g$Delta, g$omega, n0.row)
  saveRDS(results[[k]], file)
  message(sprintf('  [%3d/%3d] %-9s r = %g, Delta = %.2f, omega = %.1f  (%.1f s)',
                  k, nrow(grid), g$Test, g$r, g$Delta, g$omega,
                  as.numeric(difftime(Sys.time(), started, units = 'secs'))))
}

results <- do.call(rbind, results)

## Long format for the figures ------------------------------------------------
#
# Each panel of Figures 1, 2, 4 and 5 shows the fixed-sample design followed by
# the BSSR designs, so the fixed-sample values are added as a single extra
# category rather than repeated for every omega.

fixed.part <- unique(results[results$omega == min(OMEGAS),
                             c('Test', 'r', 'Delta', 'theta',
                               'TIE.fixed', 'power.fixed')])
names(fixed.part)[names(fixed.part) == 'TIE.fixed'] <- 'TIE'
names(fixed.part)[names(fixed.part) == 'power.fixed'] <- 'power'
fixed.part$design <- 'Fixed'

bssr.part <- results[, c('Test', 'r', 'Delta', 'theta', 'omega',
                         'TIE.BSSR', 'power.BSSR')]
names(bssr.part)[names(bssr.part) == 'TIE.BSSR'] <- 'TIE'
names(bssr.part)[names(bssr.part) == 'power.BSSR'] <- 'power'
bssr.part$design <- sprintf('%.1f', bssr.part$omega)
bssr.part$omega <- NULL

long <- rbind(fixed.part[, c('Test', 'r', 'Delta', 'theta', 'design',
                             'TIE', 'power')],
              bssr.part[, c('Test', 'r', 'Delta', 'theta', 'design',
                            'TIE', 'power')])
long$design <- factor(long$design, levels = c('Fixed', sprintf('%.1f', OMEGAS)))
long$Delta.label <- factor(sprintf('Delta = %.2f', long$Delta),
                           levels = sprintf('Delta = %.2f', sort(DELTAS, TRUE)))

## Figures --------------------------------------------------------------------

box_figure <- function(df, yvar, hline, ylab, title) {
  ggplot(df, aes(x = .data[['design']], y = .data[[yvar]])) +
    geom_boxplot(outlier.size = 0.6, linewidth = 0.3, fill = 'grey92') +
    geom_hline(yintercept = hline, linetype = 'dotted', colour = 'firebrick') +
    facet_grid(Delta.label ~ Test, scales = 'free_y') +
    labs(title = title,
         x = 'Fixed-sample design, then BSSR at information fraction omega',
         y = ylab) +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(face = 'bold'),
          panel.grid.minor = element_blank())
}

for (r in R.VALUES) {
  sub <- long[long$r == r, ]
  ggsave(
    file.path(OUT.DIR, sprintf('figure-tie-r%g.png', r)),
    box_figure(sub, 'TIE', ALPHA, 'Type I error rate',
               sprintf('Type I error rate, r = %g (Figures 1 and 2)', r)),
    width = 9, height = 10, dpi = 150
  )
  ggsave(
    file.path(OUT.DIR, sprintf('figure-power-r%g.png', r)),
    box_figure(sub, 'power', TAR.POWER, 'Power',
               sprintf('Power, r = %g (Figures 4 and 5)', r)),
    width = 9, height = 10, dpi = 150
  )
  en <- results[results$r == r, ]
  en$design <- factor(sprintf('%.1f', en$omega),
                      levels = sprintf('%.1f', OMEGAS))
  en$Delta.label <- factor(sprintf('Delta = %.2f', en$Delta),
                           levels = sprintf('Delta = %.2f', sort(DELTAS, TRUE)))
  ggsave(
    file.path(OUT.DIR, sprintf('figure-en-r%g.png', r)),
    ggplot(en, aes(x = design, y = EN.BSSR)) +
      geom_boxplot(outlier.size = 0.6, linewidth = 0.3, fill = 'grey92') +
      geom_hline(aes(yintercept = n0), linetype = 'dotted', colour = 'firebrick') +
      facet_grid(Delta.label ~ Test, scales = 'free_y') +
      labs(title = sprintf('Expected final sample size, r = %g (Figures 6 and 7)', r),
           x = 'Information fraction omega',
           y = 'Expected final sample size') +
      theme_bw(base_size = 9) +
      theme(plot.title = element_text(face = 'bold'),
            panel.grid.minor = element_blank()),
    width = 9, height = 10, dpi = 150
  )
}

## Claims made in the article -------------------------------------------------
#
# Section 4.1 and Web Appendix D report a Boschloo type I error rate of 0.025001
# for Scenario A, (r, Delta, omega, theta) = (2, 0.3, 0.9, 0.58), and a Z-pooled
# rate of 0.0252 for Scenario B, (1, 0.36, 0.8, 0.18). Section 4 reports a fixed
# sample size of n0 = 80 for Delta = 0.3 and r = 1.

lookup_tie <- function(tst, r, Delta, omega, theta) {
  n0.row <- n0.table[n0.table$Test == tst & n0.table$r == r &
                       n0.table$Delta == Delta, ]
  if (nrow(n0.row) == 0) return(NA_real_)
  run <- BinaryPowerBSSR(
    p = theta,
    Delta.A = Delta, Delta.T = 0,
    N1 = n0.row$nE0, N2 = n0.row$nC0, omega = omega, r = r,
    alpha = ALPHA, tar.power = TAR.POWER, Test = tst,
    restricted = FALSE, n.grid = N.GRID, bb.gamma = BB.GAMMA
  )
  run$power.BSSR[1]
}

claims <- data.frame(
  quantity = c('Scenario A, Boschloo TIE at (2, 0.30, 0.9, 0.58)',
               'Scenario B, Z-pool TIE at (1, 0.36, 0.8, 0.18)',
               'n0 for Delta = 0.30, r = 1, Boschloo',
               'n0 for Delta = 0.30, r = 1, Z-pool'),
  article = c(0.025001, 0.0252, 80, 80),
  reproduced = c(
    lookup_tie('Boschloo', 2, 0.30, 0.9, 0.58),
    lookup_tie('Z-pool', 1, 0.36, 0.8, 0.18),
    n0.table$n0[n0.table$Test == 'Boschloo' & n0.table$r == 1 &
                  n0.table$Delta == 0.30][1],
    n0.table$n0[n0.table$Test == 'Z-pool' & n0.table$r == 1 &
                  n0.table$Delta == 0.30][1]
  ),
  stringsAsFactors = FALSE
)

## Maximum type I error rate over the whole unit interval ---------------------
#
# The article searches the nuisance parameter over [r Delta / (1 + r),
# 1 - Delta / (1 + r)]. That restriction is what keeps both group probabilities
# inside the unit interval when the treatment effect is present, and it is
# needed for the power calculation. Under the null hypothesis the two groups
# share one probability, so the restriction is not implied by anything and the
# definition of exactness is a supremum over the whole unit interval. The block
# below repeats the type I error calculation over the full range.

message('Searching the type I error rate over the whole unit interval')

theta.full <- round(seq(0.01, 0.99, by = THETA.STEP / 2), 10)

full.range <- do.call(rbind, lapply(seq_len(nrow(grid)), function(k) {
  g <- grid[k, ]
  key <- sprintf('full_%s_r%g_d%g_w%g_g%d_bb%g',
                 g$Test, g$r, g$Delta, g$omega, N.GRID, BB.GAMMA)
  file <- file.path(OUT.DIR, 'cache', paste0(key, '.rds'))
  if (file.exists(file)) return(readRDS(file))
  n0.row <- n0.table[n0.table$Test == g$Test & n0.table$r == g$r &
                       n0.table$Delta == g$Delta, ]
  run <- BinaryPowerBSSR(
    p = theta.full,
    Delta.A = g$Delta, Delta.T = 0,
    N1 = n0.row$nE0, N2 = n0.row$nC0, omega = g$omega, r = g$r,
    alpha = ALPHA, tar.power = TAR.POWER, Test = g$Test,
    restricted = FALSE, n.grid = N.GRID, bb.gamma = BB.GAMMA
  )
  restricted.theta <- theta.grid(g$Delta, g$r)
  inside <- run$p >= min(restricted.theta) - 1e-9 &
    run$p <= max(restricted.theta) + 1e-9
  out <- data.frame(
    Test = g$Test, r = g$r, Delta = g$Delta, omega = g$omega,
    max.TIE.article.range = max(run$power.BSSR[inside]),
    max.TIE.full.range = max(run$power.BSSR),
    attained.at.theta = run$p[which.max(run$power.BSSR)],
    stringsAsFactors = FALSE
  )
  saveRDS(out, file)
  out
}))

full.range$relative.excess <- full.range$max.TIE.full.range / ALPHA - 1

## Output ---------------------------------------------------------------------

write.csv(n0.table, file.path(OUT.DIR, 'fixed-sample-size.csv'), row.names = FALSE)
write.csv(stage1.check, file.path(OUT.DIR, 'stage1-agreement.csv'), row.names = FALSE)
write.csv(results, file.path(OUT.DIR, 'operating-characteristics.csv'), row.names = FALSE)
write.csv(claims, file.path(OUT.DIR, 'article-claims.csv'), row.names = FALSE)
write.csv(full.range, file.path(OUT.DIR, 'tie-full-range.csv'), row.names = FALSE)

cat('\n== Fixed-sample size ==\n')
print(n0.table, row.names = FALSE)

cat('\n== Stage 1 sizes, settings where bbssr and the article differ ==\n')
disagree <- stage1.check[!stage1.check$agree, ]
if (nrow(disagree) == 0) {
  cat('  none, the reproduction is exact for every setting\n')
} else {
  print(disagree, row.names = FALSE)
}

cat('\n== Values reported in the article ==\n')
print(claims, row.names = FALSE)

cat('\n== Largest type I error rate, article range against the whole interval ==\n')
worst <- full.range[order(-full.range$max.TIE.full.range), ]
print(utils::head(worst, 15), row.names = FALSE)

cat(sprintf('\nSettings above the nominal level within the article range : %d of %d\n',
            sum(full.range$max.TIE.article.range > ALPHA), nrow(full.range)))
cat(sprintf('Settings above the nominal level over the whole interval  : %d of %d\n',
            sum(full.range$max.TIE.full.range > ALPHA), nrow(full.range)))
cat(sprintf('\nOutput written to %s\n', normalizePath(OUT.DIR)))
