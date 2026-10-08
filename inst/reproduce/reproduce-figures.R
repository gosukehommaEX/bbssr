# Computes with bbssr the values plotted in the published figures that the vignette
# 'Validation by reproducing published figures' redraws: Figures 1 to 4 of Friede and
# Kieser (2004), Figures 21.1 to 21.4 and 23.1 to 23.3 of Kieser (2020), Figures 1 to 3 of
# Friede, Mitchell and Mueller-Velten (2007) and Figure 2 of Boschloo (1970). The designs
# follow inst/reproduce/reproduce-published.R. The figures themselves are drawn by the
# functions in inst/reproduce/plot-figures.R.
#
# Run from the package root after devtools::load_all(); the computation for Figure 21.1
# uses internal functions of the package. The values are written as CSV files to
# inst/extdata/published-figures/, from which the vignette draws the figures. The full run
# takes about half an hour; Figure 1 of Friede et al. (2007) takes most of it.

out.dir <- file.path('inst', 'extdata', 'published-figures')
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
t.start <- Sys.time()
elapsed <- function(label) {
  message(sprintf('%s: %.1f minutes', label,
                  as.numeric(difftime(Sys.time(), t.start, units = 'mins'))))
}
save_csv <- function(x, name) {
  write.csv(x, file.path(out.dir, name), row.names = FALSE)
}

# Friede and Kieser (2004) ---------------------------------------------------------------
# Two-sided chi-squared test at level 0.05, target power 0.8, sample size from formula (1)
# (variance under the null hypothesis) with the rounding of the article, and unrestricted
# re-estimation with an upper bound of twice the required total n. For a required total n
# the assumed difference is the one for which formula (1) gives n at the overall rate pi,
# so that the re-estimated total is n p.hat (1 - p.hat) / (pi (1 - pi)) as in formula (2)
K <- (qnorm(0.975) + qnorm(0.8))^2
# Type I error rate of the fixed design with a total of n patients, allocated as
# ceiling(n / (1 + theta)) and theta times as many
fk_fixed <- function(pi, n, theta) {
  N2 <- ceiling(n / (1 + theta))
  BinaryPower(pi, pi, theta * N2, N2, 0.05, 'Chisq', alternative = 'two.sided')$Power
}
# Type I error rate of the internal pilot study design with a pilot of n.pilot patients
# and the required total n
fk_ips <- function(pi, n, n.pilot, theta) {
  n12 <- n.pilot / (1 + theta)
  N2 <- ceiling(n / (1 + theta))
  Delta <- sqrt((1 + theta)^2 / theta * K * pi * (1 - pi) / n)
  BinaryPowerBSSR(p = pi, Delta.A = Delta, Delta.T = 0, N1 = theta * N2, N2 = N2,
                  n.interim = c(theta * n12, n12), r = theta, alpha = 0.05,
                  tar.power = 0.8, Test = 'Chisq', alternative = 'two.sided',
                  ss.method = 'null.variance', rounding = 'friede-kieser',
                  N.max = max(2 * n, n.pilot))$power.BSSR
}

# Figure 1: balanced groups and the overall rates 0.1 and 0.5. For the fixed design n is
# the observed total, for the internal pilot study design the required total
fk1 <- expand.grid(n = seq(10, 200, by = 2), pi = c(0.1, 0.5))
fk1$fixed <- vapply(seq_len(nrow(fk1)), function(i) fk_fixed(fk1$pi[i], fk1$n[i], 1),
                    numeric(1))
for (n1 in seq(20, 100, by = 20)) {
  fk1[[paste0('pilot.', n1)]] <- vapply(seq_len(nrow(fk1)), function(i) {
    fk_ips(fk1$pi[i], fk1$n[i], n1, 1)
  }, numeric(1))
}
save_csv(fk1, 'friede-kieser-2004-figure1.csv')
elapsed('Friede and Kieser (2004), Figure 1')

# Figure 2: minimum, mean and maximum over the 168 combinations of the overall rate and
# the required total. The fixed design takes the observed totals of at least the pilot
fk2.grid <- expand.grid(pi = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5), n = seq(30, 300, by = 10))
fk2 <- list()
for (th in c(1, 3)) {
  fixed <- vapply(seq_len(nrow(fk2.grid)), function(i) {
    fk_fixed(fk2.grid$pi[i], fk2.grid$n[i], th)
  }, numeric(1))
  for (n1 in seq(20, 200, by = 20)) {
    ips <- vapply(seq_len(nrow(fk2.grid)), function(i) {
      fk_ips(fk2.grid$pi[i], fk2.grid$n[i], n1, th)
    }, numeric(1))
    fx <- fixed[fk2.grid$n >= n1]
    fk2[[length(fk2) + 1L]] <- data.frame(
      theta = th, pilot = n1, fixed.min = min(fx), fixed.mean = mean(fx),
      fixed.max = max(fx), ips.min = min(ips), ips.mean = mean(ips), ips.max = max(ips)
    )
  }
}
save_csv(do.call(rbind, fk2), 'friede-kieser-2004-figure2.csv')
elapsed('Friede and Kieser (2004), Figure 2')

# Figure 3: power for the alternative 0.2 as in Table I of the article. In bbssr the
# larger group is group 1, which receives the larger response probability pi1 + 0.2
fk3 <- rbind(data.frame(theta = 1, pi1 = seq(0.05, 0.40, by = 0.05)),
             data.frame(theta = 3, pi1 = seq(0.05, 0.75, by = 0.05)))
fk3$pi1 <- round(fk3$pi1, 2)
fk3[c('n', 'fixed', 'pilot.40', 'pilot.80', 'pilot.120')] <- NA_real_
for (i in seq_len(nrow(fk3))) {
  th <- fk3$theta[i]
  q1 <- fk3$pi1[i] + 0.2
  q2 <- fk3$pi1[i]
  ss <- BinarySampleSize(q1, q2, th, 0.05, 0.8, 'Chisq', alternative = 'two.sided',
                         method = 'null.variance', rounding = 'friede-kieser')
  fk3$n[i] <- ss$N
  for (n1 in c(40, 80, 120)) {
    n12 <- n1 / (1 + th)
    res <- BinaryPowerBSSR(
      p = (th * q1 + q2) / (1 + th), Delta.A = 0.2, Delta.T = 0.2, N1 = ss$N1,
      N2 = ss$N2, n.interim = c(th * n12, n12), r = th, alpha = 0.05, tar.power = 0.8,
      Test = 'Chisq', alternative = 'two.sided', ss.method = 'null.variance',
      rounding = 'friede-kieser', N.max = 2 * ss$N
    )
    fk3$fixed[i] <- res$power.TRAD
    fk3[[paste0('pilot.', n1)]][i] <- res$power.BSSR
  }
}
save_csv(fk3, 'friede-kieser-2004-figure3.csv')
elapsed('Friede and Kieser (2004), Figure 3')

# Figure 4: the depression trial of Section 5. Left, the fixed design with 122 patients
# per group; right, the internal pilot study design with the alternative 0.15 and a pilot
# of 120 patients against the required total n, converted to the overall rate with
# formula (1)
fk4a <- data.frame(pi = seq(0.01, 0.50, by = 0.001))
fk4a$level <- BinaryPower(fk4a$pi, fk4a$pi, 122, 122, 0.05, 'Chisq',
                          alternative = 'two.sided')$Power
save_csv(fk4a, 'friede-kieser-2004-figure4-fixed.csv')
fk4b <- data.frame(n = seq(20, 348, by = 2))
fk4b$pi <- (1 - sqrt(1 - 4 * fk4b$n * 0.15^2 / (4 * K))) / 2
fk4b$level <- vapply(seq_len(nrow(fk4b)), function(i) {
  BinaryPowerBSSR(p = fk4b$pi[i], Delta.A = 0.15, Delta.T = 0, N1 = 122, N2 = 122,
                  n.interim = c(60, 60), r = 1, alpha = 0.05, tar.power = 0.8,
                  Test = 'Chisq', alternative = 'two.sided', ss.method = 'null.variance',
                  rounding = 'friede-kieser', N.max = max(2 * fk4b$n[i], 120))$power.BSSR
}, numeric(1))
save_csv(fk4b, 'friede-kieser-2004-figure4-recalculation.csv')
elapsed('Friede and Kieser (2004), Figure 4')

# Kieser (2020) ------------------------------------------------------------------------
# Normal approximation test (equivalent to the chi-squared test) at the one-sided level
# 0.025, target power 0.8, sample size from formula (21.3) and unrestricted
# re-estimation. Group 1 of bbssr is the experimental group E, so the allocation ratio
# r = nE / nC of the book is the ratio r of bbssr

# Type I error rate of the recalculation design when the trial stops after the pilot for
# the interim outcomes whose recovered rate of a group lies outside the unit interval. The
# final sample sizes of bbssr are kept for all other outcomes; bbssr itself truncates the
# Bernoulli variance of such a group at zero and re-estimates the sample size
kieser_pilot_stop <- function(p, Delta, r, N1, N2, n.interim) {
  map <- bssr_map(Delta.A = Delta, N1 = N1, N2 = N2, omega = NULL, n.interim = n.interim,
                  r = r, alpha = 0.025, tar.power = 0.8, Test = 'Chisq',
                  restricted = FALSE, alternative = 'greater', tsmethod = 'minlike',
                  n.grid = 100, bb.gamma = 0, effect = 'RD', ss.method = 'standard',
                  ss.Test = 'Chisq', ss.alpha = 0.025, rounding = 'group', N.min = NULL,
                  N.max = NULL, ref.pvalue = FALSE, margin = 0)
  sp <- split_pooled(map$hat.p, Delta, r, 'RD')
  outside <- sp$p1 > 1 + 1e-9 | sp$p2 < -1e-9
  map$N1[outside] <- attr(map, 'n11')
  map$N2[outside] <- attr(map, 'n12')
  setup <- bssr_setup(map)
  rr.list <- lapply(seq_along(setup$N1), function(k) {
    get_rr(setup$N1[k], setup$N2[k], 0.025, 'Chisq', 'greater', 'minlike', 100, 0, FALSE,
           0)
  })
  bssr_reject(setup, rr.list, p, p)
}

# Figure 21.1: actual levels at the assumed overall rate for the alternatives 0.15 to 0.30
# (pE > pC), the assumed overall rates 0.30 to 0.50 and pilots of 25, 50 and 75 per cent
# of the fixed sample size, with ceiling(fraction nC) patients in group C and r times as
# many in group E (79 + 79 in Example 21.1)
k1 <- expand.grid(pA = round(seq(0.30, 0.50, by = 0.01), 2),
                  Delta = c(0.15, 0.20, 0.25, 0.30), fraction = c(0.25, 0.5, 0.75),
                  r = c(1, 3))
k1[c('nE', 'nC', 'n1E', 'n1C', 'fixed', 'ips', 'ips.stop')] <- NA_real_
for (i in seq_len(nrow(k1))) {
  ri <- k1$r[i]
  pE <- k1$pA[i] + k1$Delta[i] / (1 + ri)
  pC <- k1$pA[i] - ri * k1$Delta[i] / (1 + ri)
  ss <- BinarySampleSize(pE, pC, ri, 0.025, 0.8, 'Chisq', method = 'standard')
  n1C <- ceiling(k1$fraction[i] * ss$N2)
  res <- BinaryPowerBSSR(p = k1$pA[i], Delta.A = k1$Delta[i], Delta.T = 0, N1 = ss$N1,
                         N2 = ss$N2, n.interim = c(ri * n1C, n1C), r = ri, alpha = 0.025,
                         tar.power = 0.8, Test = 'Chisq', ss.method = 'standard')
  stop.level <- kieser_pilot_stop(k1$pA[i], k1$Delta[i], ri, ss$N1, ss$N2,
                                  c(ri * n1C, n1C))
  k1[i, c('nE', 'nC', 'n1E', 'n1C', 'fixed', 'ips', 'ips.stop')] <-
    c(ss$N1, ss$N2, ri * n1C, n1C, res$power.TRAD, res$power.BSSR, stop.level)
}
save_csv(k1, 'kieser-2020-figure21-1.csv')
elapsed('Kieser (2020), Figure 21.1')

# Figure 21.2: Example 21.1 (Delta = -0.15, n0 = 2 x 157, pilot of 79 + 79) and the same
# design with Delta = -0.30 (n0 = 2 x 39, pilot of 20 + 20) on the grid of the book. The
# figure shows the overall rates for which the assumed rates of both groups lie in the
# unit interval, [0.075, 0.925] and [0.15, 0.85]
k.args <- list(Delta.A = -0.15, N1 = 157, N2 = 157, n.interim = c(79, 79), r = 1,
               alpha = 0.025, tar.power = 0.8, Test = 'Chisq', alternative = 'less',
               ss.method = 'standard', theta = seq(0.01, 0.99, by = 0.005),
               maximize = 'grid')
k2 <- do.call(rbind, lapply(list(k.args, modifyList(k.args, list(
  Delta.A = -0.30, N1 = 39, N2 = 39, n.interim = c(20, 20)
))), function(a) {
  x <- do.call(BinaryTypeIErrorBSSR, a)
  data.frame(Delta = a$Delta.A, p = x$theta, fixed = x$TIE.TRAD,
             recalculation = x$TIE.BSSR)
}))
k2$shown <- k2$p >= abs(k2$Delta) / 2 - 1e-9 & k2$p <= 1 - abs(k2$Delta) / 2 + 1e-9
save_csv(k2, 'kieser-2020-figure21-2.csv')
elapsed('Kieser (2020), Figure 21.2')

# Figure 21.3: power for the rates pE = 0.27 and pC = 0.42 of Example 21.1 with r = 1 and
# r = 3. For r = 3 the book does not state the sample sizes; formula (21.3) gives them for
# the same rates (overall rate 0.3075), and the pilot is half of each group, rounded as in
# Figure 21.1
k3 <- do.call(rbind, lapply(c(1, 3), function(ri) {
  ss <- BinarySampleSize(0.27, 0.42, ri, 0.025, 0.8, 'Chisq', alternative = 'less',
                         method = 'standard')
  n1C <- ceiling(ss$N2 / 2)
  # Overall rates for which both true rates lie in the unit interval
  p <- seq(0.15 / (1 + ri), 1 - ri * 0.15 / (1 + ri), by = 0.005)
  res <- BinaryPowerBSSR(p = p, Delta.A = -0.15, Delta.T = -0.15, N1 = ss$N1, N2 = ss$N2,
                         n.interim = c(ri * n1C, n1C), r = ri, alpha = 0.025,
                         tar.power = 0.8, Test = 'Chisq', alternative = 'less',
                         ss.method = 'standard')
  data.frame(r = ri, nE = ss$N1, nC = ss$N2, n1E = ri * n1C, n1C = n1C, p = res$p,
             fixed = res$power.TRAD, recalculation = res$power.BSSR)
}))
save_csv(k3, 'kieser-2020-figure21-3.csv')
elapsed('Kieser (2020), Figure 21.3')

# Box plots of the exact distribution of the recalculated total sample size for each
# overall rate of res, a result of BinaryPowerBSSR(). The quartiles are the smallest
# totals whose cumulative probability reaches them, as in summary(); the whiskers reach
# the most extreme totals within 1.5 interquartile ranges of the box, and the totals
# beyond them with a probability of at least 1e-4 are kept for plotting
kieser_boxes <- function(res) {
  k4.box <- summary(res)[c('p', 'E.N', 'N.q25', 'N.q50', 'N.q75')]
  k4.box[c('whisker.low', 'whisker.high')] <- NA_real_
  k4.dist <- attr(res, 'N.dist')
  k4.out <- list(data.frame(p = numeric(0), N = numeric(0), prob = numeric(0)))
  for (k in seq_len(nrow(k4.box))) {
    d <- k4.dist[k4.dist$scenario == k, ]
    prob <- tapply(d$prob, d$N, sum)
    N <- as.numeric(names(prob))
    iqr <- k4.box$N.q75[k] - k4.box$N.q25[k]
    inside <- N >= k4.box$N.q25[k] - 1.5 * iqr & N <= k4.box$N.q75[k] + 1.5 * iqr
    k4.box$whisker.low[k] <- min(N[inside])
    k4.box$whisker.high[k] <- max(N[inside])
    far <- !inside & prob >= 1e-4
    if (any(far)) {
      k4.out[[length(k4.out) + 1L]] <- data.frame(
        p = k4.box$p[k], N = N[far], prob = as.numeric(prob[far])
      )
    }
  }
  list(box = k4.box, outside = do.call(rbind, k4.out))
}

# Figure 21.4: distribution of the recalculated total sample size of Example 21.1, with
# the total rounded up as in inst/reproduce/reproduce-published.R
k4 <- BinaryPowerBSSR(p = seq(0.10, 0.90, by = 0.05), Delta.A = -0.15, Delta.T = -0.15,
                      N1 = 157, N2 = 157, n.interim = c(79, 79), r = 1, alpha = 0.025,
                      tar.power = 0.8, Test = 'Chisq', alternative = 'less',
                      ss.method = 'standard', rounding = 'total')
k4.boxes <- kieser_boxes(k4)
save_csv(k4.boxes$box, 'kieser-2020-figure21-4.csv')
save_csv(k4.boxes$outside, 'kieser-2020-figure21-4-outside.csv')
elapsed('Kieser (2020), Figure 21.4')

# Example 23.1 (FreezeAF trial): Farrington-Manning test of pE - pC <= -0.15 at the
# one-sided level 0.025, the equal assumed rates 0.78, power 0.8, n0 = 2 x 122 from
# formula (23.2) with each group rounded up, a pilot of 50 + 50 patients and unrestricted
# recalculation, as in inst/reproduce/reproduce-published.R
k23.args <- list(Delta.A = 0, N1 = 122, N2 = 122, n.interim = c(50, 50), r = 1,
                 alpha = 0.025, tar.power = 0.8, Test = 'Farrington-Manning',
                 ss.method = 'standard', margin = 0.15)
k23.p <- round(seq(0.075, 0.925, by = 0.005), 3)

# Figure 23.1: actual level on the boundary of the null hypothesis, pE = p - 0.075 and
# pC = p + 0.075, against the overall rate p over [0.075, 0.925]
k23.1 <- do.call(BinaryTypeIErrorBSSR,
                 c(k23.args, list(theta = k23.p, maximize = 'grid')))
save_csv(data.frame(p = k23.1$theta, fixed = k23.1$TIE.TRAD,
                    recalculation = k23.1$TIE.BSSR), 'kieser-2020-figure23-1.csv')
elapsed('Kieser (2020), Figure 23.1')

# Figure 23.2: power for the equal true rates pE = pC = p
k23.2 <- do.call(BinaryPowerBSSR, c(k23.args, list(p = k23.p, Delta.T = 0)))
save_csv(data.frame(p = k23.2$p, fixed = k23.2$power.TRAD,
                    recalculation = k23.2$power.BSSR), 'kieser-2020-figure23-2.csv')
elapsed('Kieser (2020), Figure 23.2')

# Figure 23.3: distribution of the recalculated total sample size for the equal true rates
# pE = pC = p
k23.3 <- kieser_boxes(do.call(BinaryPowerBSSR, c(k23.args, list(
  p = round(seq(0.10, 0.90, by = 0.05), 2), Delta.T = 0
))))
save_csv(k23.3$box, 'kieser-2020-figure23-3.csv')
save_csv(k23.3$outside, 'kieser-2020-figure23-3-outside.csv')
elapsed('Kieser (2020), Figure 23.3')

# Boschloo (1970), Figure 2 --------------------------------------------------------------
# Fisher's test of p1 <= p2 against p1 > p2 at the conditional level 0.05 with m = 15 and
# n = 10 (N1 = m, N2 = n). For each total number of successes r, alpha_r is the
# conditional size given r, and f_r(p) = alpha_r p(r | H0) is the contribution of r to the
# unconditional level alpha(p). alpha_gamma(p) is the level of Fisher's test at the
# raised conditional level gamma = 0.09, which is the Boschloo test at level 0.05
bo.cells <- matrix(as.numeric(BinaryRR(15, 10, 0.05, 'Fisher')), 16, 11)
bo.a <- row(bo.cells) - 1
bo.b <- col(bo.cells) - 1
bo.alpha.r <- vapply(0:25, function(r) {
  k <- bo.a + bo.b == r
  sum(bo.cells[k] * dhyper(bo.a[k], 15, 10, r))
}, numeric(1))
bo.p <- seq(0, 1, by = 0.002)
bo.f <- vapply(0:25, function(r) bo.alpha.r[r + 1] * dbinom(r, 25, bo.p),
               numeric(length(bo.p)))
bo <- data.frame(p = bo.p, alpha = rowSums(bo.f),
                 alpha.gamma = BinaryPower(bo.p, bo.p, 15, 10, 0.09, 'Fisher')$Power,
                 f6 = bo.f[, 7], f10 = bo.f[, 11], f17 = bo.f[, 18], f21 = bo.f[, 22])
# The sum of the f_r is the level that BinaryPower() computes from the rejection region,
# and Fisher's test at the raised level is the Boschloo test
bo.fisher <- BinaryPower(bo.p, bo.p, 15, 10, 0.05, 'Fisher')$Power
bo.boschloo <- BinaryPower(bo.p, bo.p, 15, 10, 0.05, 'Boschloo')$Power
stopifnot(max(abs(bo$alpha - bo.fisher)) < 1e-12,
          max(abs(bo$alpha.gamma - bo.boschloo)) < 1e-12)
save_csv(bo, 'boschloo-1970-figure2.csv')
save_csv(data.frame(r = 0:25, alpha.r = bo.alpha.r, p = (0:25) / 25,
                    f.max = bo.alpha.r * dbinom(0:25, 25, (0:25) / 25)),
         'boschloo-1970-figure2-maxima.csv')
elapsed('Boschloo (1970), Figure 2')

# Friede, Mitchell and Mueller-Velten (2007) ---------------------------------------------
# One-sided level 0.025, margin 0.1, target power 0.8 and assumed difference 0. Group 1 is
# the experimental group and the article writes q = n2 / n1, so r = 1 / q. The fixed
# sample sizes come from formula (1) (Blackwelder) or (2) (Farrington-Manning) with the
# two groups rounded up separately, the internal pilot study has half of each group
# rounded to the nearest integer, and the re-estimated groups are rounded to the nearest
# integer. The article writes delta1 = p2 - p1, which is -Delta.T in bbssr
fmm_method <- function(Test) {
  if (Test == 'Blackwelder') 'alternative.variance' else 'standard'
}
fmm_design <- function(pa, q, Test, p, delta1) {
  ss <- BinarySampleSize(pa, pa, 1 / q, 0.025, 0.8, Test, method = fmm_method(Test),
                         rounding = 'friede-kieser', margin = 0.1)
  BinaryPowerBSSR(p = p, Delta.A = 0, Delta.T = -delta1, N1 = ss$N1, N2 = ss$N2,
                  n.interim = floor(c(ss$N1, ss$N2) / 2 + 0.5), r = 1 / q, alpha = 0.025,
                  tar.power = 0.8, Test = Test, ss.method = fmm_method(Test),
                  rounding = 'nearest', margin = 0.1)
}
shift <- c(-0.2, -0.1, 0, 0.1, 0.2)

# Figure 2: power of the Farrington-Manning test with the correctly specified difference
# 0, for pa = 0.5 and 0.7, q = 1 and 1/3 and the true overall rates pa - 0.2 to pa + 0.2
m2 <- do.call(rbind, lapply(c(1, 1 / 3), function(q) {
  do.call(rbind, lapply(c(0.5, 0.7), function(pa) {
    res <- fmm_design(pa, q, 'Farrington-Manning', pa + shift, 0)
    data.frame(q = q, pa = pa, shift = round(res$p - pa, 2), fixed = res$power.TRAD,
               reestimation = res$power.BSSR)
  }))
}))
save_csv(m2, 'friede-2007-figure2.csv')
elapsed('Friede et al. (2007), Figure 2')

# Figure 3: power of the Farrington-Manning test with q = 1 when the true difference
# delta1 differs from the assumed difference 0 by -0.02 to 0.02
m3 <- do.call(rbind, lapply(c(0.5, 0.7), function(pa) {
  do.call(rbind, lapply(c(-0.02, -0.01, 0, 0.01, 0.02), function(d1) {
    res <- fmm_design(pa, 1, 'Farrington-Manning', pa + shift, d1)
    data.frame(pa = pa, delta1 = d1, shift = round(res$p - pa, 2), fixed = res$power.TRAD,
               reestimation = res$power.BSSR)
  }))
}))
save_csv(m3, 'friede-2007-figure3.csv')
elapsed('Friede et al. (2007), Figure 3')

# Figure 1: type I error rates on the boundary of the null hypothesis (delta1 = 0.1) at
# the assumed overall rate, for q = 1/3, 1/2 and 1 and the assumed overall rates 0.30 to
# 0.70
m1 <- expand.grid(pa = round(seq(0.30, 0.70, by = 0.01), 2), q = c(1 / 3, 1 / 2, 1),
                  Test = c('Blackwelder', 'Farrington-Manning'), stringsAsFactors = FALSE)
m1[c('N1', 'N2', 'fixed', 'reestimation')] <- NA_real_
for (i in seq_len(nrow(m1))) {
  res <- fmm_design(m1$pa[i], m1$q[i], m1$Test[i], m1$pa[i], 0.1)
  m1[i, c('N1', 'N2', 'fixed', 'reestimation')] <-
    c(attr(res, 'N1'), attr(res, 'N2'), res$power.TRAD, res$power.BSSR)
  if (i %% 41 == 0) elapsed(sprintf('Friede et al. (2007), Figure 1, %d of %d designs', i,
                                    nrow(m1)))
}
save_csv(m1, 'friede-2007-figure1.csv')
elapsed('Friede et al. (2007), Figure 1')
