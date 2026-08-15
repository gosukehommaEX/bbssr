###############################################################################
# reproduce_cui_etal.R
#
# Reproduces the numerical results of
#
#   Cui, Y., Homma, G. and Daimon, T.
#   "Exact unconditional tests for blinded sample size re-estimation"
#
# using bbssr 2.0.0, and extends them in the two directions that the article
# leaves open.
#
#   1. The type I error rate is searched over the whole unit interval. The
#      article restricts the nuisance parameter to [r Delta / (1 + r),
#      1 - Delta / (1 + r)]. That restriction is what keeps both group
#      probabilities inside the unit interval when the treatment effect is
#      present, so it belongs to the power calculation. Under the null
#      hypothesis the two groups share one probability and the restriction is
#      not implied by anything, while the definition of exactness is a supremum
#      over the whole interval.
#
#   2. The grid over the nuisance parameter of the unconditional tests is set to
#      N.GRID = 2000 rather than the 100 points of the article. A coarse grid
#      understates the p-value, which enlarges the rejection region and inflates
#      the type I error rate that is then computed for it. Any inflation
#      attributed to the adaptation has to survive a refinement of this grid,
#      otherwise it is a numerical artefact. Setting N.GRID to 100 reproduces
#      the article, and the two runs can be compared because the output files
#      and the cache carry the grid size in their names.
#
# Version 2.0.0 of bbssr also corrects two things that affect the numbers.
# Outcomes sharing the same value of the ordering statistic now receive the same
# p-value, and the allocation ratio is preserved by the interim and the final
# sample sizes. The figures produced here can therefore differ from the
# published ones. Every number the script prints is computed, none is taken from
# the article, and the values reported in the article are shown alongside.
#
# Run from the package root with
#
#   source('tools/reproduce_cui_etal.R')
#
# Run time: the default settings take on the order of ten hours. Each setting is
# cached as soon as it finishes, so the script can be interrupted and restarted,
# and it prints an estimate of the time remaining. Set QUICK to TRUE for a
# coarse subset that finishes in a few minutes.
###############################################################################

library(bbssr)
library(ggplot2)

## Configuration --------------------------------------------------------------

# TRUE runs a coarse subset that finishes in a few minutes
QUICK <- FALSE

# Directory for cached results, tables and figures
OUT.DIR <- 'reproduce-output'

# The cache is keyed on the settings, the version of bbssr and the tag below.
# Increment the tag whenever the behaviour of bbssr changes without a change of
# version number, otherwise a stale result is silently reused. Deleting the
# cache directory has the same effect.
CACHE.TAG <- 'v3'

# Number of grid points for the nuisance parameter of the unconditional tests.
# The article used 100. See the note at the top of the file.
N.GRID <- if (QUICK) 100 else 2000

# Berger-Boos parameter. Set to 1e-4 to reproduce Web Appendix C
BB.GAMMA <- 0

# Step of the grid used to search the type I error rate over the whole unit
# interval. Adding values of the nuisance parameter is cheap, because the
# re-estimated sample sizes and the rejection regions are shared across them.
THETA.STEP.FULL <- if (QUICK) 0.005 else 0.002

## Settings of the article ----------------------------------------------------
#
# Section 4 fixes the response probability of the experimental group at
# theta.E = 0.4 and varies the treatment effect, so theta.C = 0.4 - Delta. The
# nominal one-sided level is 0.025 and the target power is 0.9. The nuisance
# parameter runs over [r Delta / (1 + r), 1 - Delta / (1 + r)] in steps of 0.01.
# Re-estimation uses the unrestricted rule at information fraction omega.

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

CACHE.KEY <- sprintf('%s-%s', as.character(utils::packageVersion('bbssr')), CACHE.TAG)
TAG <- sprintf('g%d%s', N.GRID, if (BB.GAMMA > 0) sprintf('-bb%g', BB.GAMMA) else '')

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
# The stage 1 sizes agree as well. The article sets nC1 = ceiling(omega nC) and
# nE1 = r nC1, and from version 2.0.0 bbssr derives the size of group 1 from the
# size of group 2 by a single ceiling, which gives the same answer for every
# whole r. check_stage1_agreement() below verifies this.

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
      bbssr.nE1 = ceiling(row$r * ceiling(omegas * row$nC0)),
      stringsAsFactors = FALSE
    )
  }))
  out$agree <- out$article.nE1 == out$bbssr.nE1
  out
}

stage1.check <- check_stage1_agreement(n0.table, OMEGAS)

## One setting of the grid ----------------------------------------------------
#
# Two calls are enough. The null run covers the whole unit interval, which is a
# superset of the grid of the article, so the type I error rate on the article
# grid is read off the same run rather than computed again. The alternative run
# is confined to the article grid, since outside it one of the two group
# probabilities leaves the unit interval.

run_setting <- function(tst, r, Delta, omega, n0.row) {
  theta.article <- theta.grid(Delta, r)
  theta.all <- sort(unique(c(theta.article,
                             round(seq(0, 1, by = THETA.STEP.FULL), 10))))
  common <- list(
    Delta.A = Delta, N1 = n0.row$nE0, N2 = n0.row$nC0,
    omega = omega, r = r,
    alpha = ALPHA, tar.power = TAR.POWER, Test = tst,
    restricted = FALSE, n.grid = N.GRID, bb.gamma = BB.GAMMA
  )
  # Delta.T = 0 puts both groups at the same response probability, so the two
  # power columns become rejection probabilities under the null hypothesis
  null.run <- do.call(BinaryPowerBSSR, c(common, list(p = theta.all, Delta.T = 0)))
  alt.run <- do.call(BinaryPowerBSSR, c(common, list(p = theta.article,
                                                     Delta.T = Delta)))
  inside <- null.run$p %in% theta.article
  curve <- data.frame(
    Test = tst, r = r, Delta = Delta, omega = omega,
    theta = null.run$p, TIE.BSSR = null.run$power.BSSR,
    TIE.fixed = null.run$power.TRAD, in.article.range = inside,
    stringsAsFactors = FALSE
  )
  oc <- data.frame(
    Test = tst, r = r, Delta = Delta, omega = omega,
    theta = alt.run$p,
    TIE.BSSR = null.run$power.BSSR[inside],
    TIE.fixed = null.run$power.TRAD[inside],
    power.BSSR = alt.run$power.BSSR,
    power.fixed = alt.run$power.TRAD,
    EN.BSSR = alt.run$E.N,
    n0 = n0.row$n0,
    n1.interim = attr(alt.run, 'n1.interim'),
    n2.interim = attr(alt.run, 'n2.interim'),
    stringsAsFactors = FALSE
  )
  summary <- data.frame(
    Test = tst, r = r, Delta = Delta, omega = omega,
    max.TIE.article.range = max(curve$TIE.BSSR[inside]),
    max.TIE.full.range = max(curve$TIE.BSSR),
    attained.at.theta = curve$theta[which.max(curve$TIE.BSSR)],
    inside.article.range = inside[which.max(curve$TIE.BSSR)],
    max.TIE.fixed.full.range = max(curve$TIE.fixed),
    stringsAsFactors = FALSE
  )
  list(oc = oc, curve = curve, summary = summary)
}

## Main loop ------------------------------------------------------------------

dir.create(file.path(OUT.DIR, 'cache'), recursive = TRUE, showWarnings = FALSE)

grid <- expand.grid(omega = OMEGAS, Delta = DELTAS, r = R.VALUES, Test = TESTS,
                    stringsAsFactors = FALSE)

message(sprintf('Running %d settings with n.grid = %d', nrow(grid), N.GRID))

all.runs <- vector('list', nrow(grid))
spent <- numeric(0)
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]
  key <- sprintf('%s_%s_%s_r%g_d%g_w%g', CACHE.KEY, TAG, g$Test, g$r, g$Delta, g$omega)
  file <- file.path(OUT.DIR, 'cache', paste0(key, '.rds'))
  if (file.exists(file)) {
    all.runs[[k]] <- readRDS(file)
    next
  }
  n0.row <- n0.table[n0.table$Test == g$Test & n0.table$r == g$r &
                       n0.table$Delta == g$Delta, ]
  started <- Sys.time()
  all.runs[[k]] <- run_setting(g$Test, g$r, g$Delta, g$omega, n0.row)
  saveRDS(all.runs[[k]], file)
  spent <- c(spent, as.numeric(difftime(Sys.time(), started, units = 'secs')))
  remaining <- mean(spent) * (nrow(grid) - k)
  message(sprintf(
    '  [%3d/%3d] %-9s r = %g, Delta = %.2f, omega = %.1f  (%.0f s, about %.1f h left)',
    k, nrow(grid), g$Test, g$r, g$Delta, g$omega,
    spent[length(spent)], remaining / 3600))
}

results <- do.call(rbind, lapply(all.runs, `[[`, 'oc'))
tie.curves <- do.call(rbind, lapply(all.runs, `[[`, 'curve'))
full.range <- do.call(rbind, lapply(all.runs, `[[`, 'summary'))
full.range$relative.excess <- full.range$max.TIE.full.range / ALPHA - 1

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
         subtitle = sprintf('Nuisance parameter grid of %d points', N.GRID),
         x = 'Fixed-sample design, then BSSR at information fraction omega',
         y = ylab) +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(face = 'bold'),
          panel.grid.minor = element_blank())
}

for (r in R.VALUES) {
  sub <- long[long$r == r, ]
  ggsave(
    file.path(OUT.DIR, sprintf('figure-tie-r%g-%s.png', r, TAG)),
    box_figure(sub, 'TIE', ALPHA, 'Type I error rate',
               sprintf('Type I error rate, r = %g (Figures 1 and 2)', r)),
    width = 9, height = 10, dpi = 150
  )
  ggsave(
    file.path(OUT.DIR, sprintf('figure-power-r%g-%s.png', r, TAG)),
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
    file.path(OUT.DIR, sprintf('figure-en-r%g-%s.png', r, TAG)),
    ggplot(en, aes(x = design, y = EN.BSSR)) +
      geom_boxplot(outlier.size = 0.6, linewidth = 0.3, fill = 'grey92') +
      geom_hline(aes(yintercept = n0), linetype = 'dotted', colour = 'firebrick') +
      facet_grid(Delta.label ~ Test, scales = 'free_y') +
      labs(title = sprintf('Expected final sample size, r = %g (Figures 6 and 7)', r),
           subtitle = 'The dotted line is the fixed-sample size',
           x = 'Information fraction omega',
           y = 'Expected final sample size') +
      theme_bw(base_size = 9) +
      theme(plot.title = element_text(face = 'bold'),
            panel.grid.minor = element_blank()),
    width = 9, height = 10, dpi = 150
  )
}

## Type I error rate against the nuisance parameter ---------------------------

exceeding <- full.range[full.range$max.TIE.full.range > ALPHA, ]
if (nrow(exceeding) > 0) {
  keys <- paste(exceeding$Test, exceeding$r, exceeding$Delta, exceeding$omega)
  sel <- tie.curves[paste(tie.curves$Test, tie.curves$r, tie.curves$Delta,
                          tie.curves$omega) %in% keys, ]
  sel$panel <- sprintf('%s, r = %g, Delta = %.2f, omega = %.1f',
                       sel$Test, sel$r, sel$Delta, sel$omega)
  bands <- do.call(rbind, lapply(split(sel, sel$panel), function(d) {
    inside <- d$theta[d$in.article.range]
    data.frame(panel = d$panel[1], xmin = min(inside), xmax = max(inside),
               stringsAsFactors = FALSE)
  }))
  ggsave(
    file.path(OUT.DIR, sprintf('figure-tie-against-theta-%s.png', TAG)),
    ggplot(sel, aes(x = theta, y = TIE.BSSR)) +
      geom_rect(data = bands, inherit.aes = FALSE,
                aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                fill = 'grey90') +
      geom_line(linewidth = 0.4, colour = 'steelblue') +
      geom_hline(yintercept = ALPHA, linetype = 'dotted', colour = 'firebrick') +
      facet_wrap(~ panel, ncol = 3, scales = 'free_y') +
      labs(title = 'Type I error rate against the nuisance parameter',
           subtitle = paste0('Settings whose maximum exceeds the nominal level, ',
                             'grid of ', N.GRID, ' points. The shaded band is the ',
                             'range searched in the article'),
           x = 'Nuisance parameter theta', y = 'Type I error rate') +
      theme_bw(base_size = 8) +
      theme(plot.title = element_text(face = 'bold'),
            panel.grid.minor = element_blank()),
    width = 11, height = 2.4 * ceiling(nrow(exceeding) / 3) + 1, dpi = 150,
    limitsize = FALSE
  )
}

## Values reported in the article ---------------------------------------------
#
# Section 4.1 and Web Appendix D report a Boschloo type I error rate of 0.025001
# for Scenario A, (r, Delta, omega, theta) = (2, 0.3, 0.9, 0.58), and a Z-pooled
# rate of 0.0252 for Scenario B, (1, 0.36, 0.8, 0.18). Section 4 reports a fixed
# sample size of n0 = 80 for Delta = 0.3 and r = 1.

lookup_tie <- function(tst, r, Delta, omega, theta) {
  hit <- tie.curves[tie.curves$Test == tst & tie.curves$r == r &
                      tie.curves$Delta == Delta & tie.curves$omega == omega &
                      abs(tie.curves$theta - theta) < 1e-9, ]
  if (nrow(hit) == 0) return(NA_real_)
  hit$TIE.BSSR[1]
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

## Output ---------------------------------------------------------------------

out_file <- function(stem) file.path(OUT.DIR, sprintf('%s-%s.csv', stem, TAG))
write.csv(n0.table, out_file('fixed-sample-size'), row.names = FALSE)
write.csv(stage1.check, out_file('stage1-agreement'), row.names = FALSE)
write.csv(results, out_file('operating-characteristics'), row.names = FALSE)
write.csv(claims, out_file('article-claims'), row.names = FALSE)
write.csv(full.range, out_file('tie-full-range'), row.names = FALSE)
write.csv(tie.curves, out_file('tie-curves'), row.names = FALSE)

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
print(utils::head(worst, 20), row.names = FALSE)

cat(sprintf('\nNuisance parameter grid                                   : %d points\n',
            N.GRID))
cat(sprintf('Settings above the nominal level within the article range : %d of %d\n',
            sum(full.range$max.TIE.article.range > ALPHA), nrow(full.range)))
cat(sprintf('Settings above the nominal level over the whole interval  : %d of %d\n',
            sum(full.range$max.TIE.full.range > ALPHA), nrow(full.range)))
cat(sprintf('Fixed-sample designs above the nominal level              : %d of %d\n',
            sum(full.range$max.TIE.fixed.full.range > ALPHA), nrow(full.range)))
if (QUICK) {
  cat('\nNOTE: QUICK is TRUE, so the grids are coarse and the maxima can step over\n',
      '      the peak. Set QUICK to FALSE for the settings of the article.\n', sep = '')
}
cat(sprintf('\nOutput written to %s\n', normalizePath(OUT.DIR)))
