# Reproduces the published numerical results of Friede and Kieser (2004, Table I and
# Section 5), Kieser (2020, Example 21.1 and the values stated for Figures 21.1 and 21.2),
# Farrington and Manning (1990, Tables I and II and the first example), Blackwelder (1982,
# Table 3 and the examples), Friede, Mitchell and Mueller-Velten (2007, Tables 2 and 3 and
# Sections 5 and 6), Boschloo (1970, Sections 2, 3, 4 and 6), Mehrotra, Chan and Berger
# (2003, Section 3.1 and Tables 1 and 3), Berger and Boos (1994, Example 2) and Fay and
# Hunsberger (2021, Section 8 and Table 1) with bbssr.
# Run from the package root after devtools::load_all(); the levels of Figure 21.1 of
# Kieser (2020) under the rule that stops after the pilot use internal functions of the
# package. The results are written to reproduce-output/: published-comparison.csv lists
# every published value next to the recomputed value with a verdict (PASS, EXPLAINED,
# FAIL or INFO), and summary.md is generated from it. The full run takes several minutes.

out.dir <- 'reproduce-output'
dir.create(out.dir, showWarnings = FALSE)
t.start <- Sys.time()
cmp <- list()
add <- function(source, item, published, recomputed, digits, note = '', rule = 'round') {
  cmp[[length(cmp) + 1L]] <<- data.frame(source = source, item = item,
                                         published = published, recomputed = recomputed,
                                         digits = digits, rule = rule, note = note,
                                         stringsAsFactors = FALSE)
}
# Discrepancies that are understood. Each entry lists the items it covers, the largest
# absolute difference from the published value that its reason accounts for, and the
# reason. An item that does not pass is EXPLAINED only when an entry covers it and its
# difference is within the tolerance of that entry; otherwise it is a FAIL, so an entry
# never hides an unexpected difference
explained <- list()
explain <- function(items, tol, reason) {
  explained[[length(explained) + 1L]] <<- list(items = items, tol = tol, reason = reason)
}

# Friede and Kieser (2004), Table I -----------------------------------------------------
# Two-sided chi-squared test at level 0.05, target power 0.8, alternative 0.2, sample size
# from formula (1) (variance under the null hypothesis) with the rounding of the article,
# unrestricted re-estimation with an upper bound of twice the initial total. In bbssr the
# larger group is group 1, which receives the larger response probability p1 + 0.2
fk <- data.frame(
  theta = rep(c(1, 3), times = c(9, 15)),
  p1 = c(rep(c(0.05, 0.2, 0.4), each = 3), rep(c(0.05, 0.2, 0.4, 0.6, 0.75), each = 3)),
  n.pilot = rep(c(40, 80, 120), 8),
  n.pub = c(rep(c(102, 166, 198), each = 3), rep(c(168, 239, 260, 198, 95), each = 3)),
  EN.pub = c(99.1, 101.8, 120.9, 162.0, 164.0, 164.6, 192.5, 195.0, 195.9,
             164.7, 166.6, 167.5, 233.5, 236.5, 237.4, 253.8, 257.1, 258.1,
             192.7, 195.1, 195.9, 94.7, 99.5, 121.3),
  pw.pub = c(0.779, 0.818, 0.898, 0.794, 0.802, 0.806, 0.802, 0.810, 0.813,
             0.846, 0.864, 0.870, 0.818, 0.825, 0.827, 0.796, 0.803, 0.803,
             0.765, 0.775, 0.779, 0.694, 0.751, 0.846)
)
fk$n <- NA
fk$EN <- NA
fk$power <- NA
for (i in seq_len(nrow(fk))) {
  th <- fk$theta[i]
  q1 <- fk$p1[i] + 0.2
  q2 <- fk$p1[i]
  ss <- BinarySampleSize(q1, q2, th, 0.05, 0.8, 'Chisq', alternative = 'two.sided',
                         method = 'null.variance', rounding = 'friede-kieser')
  n12 <- fk$n.pilot[i] / (1 + th)
  res <- BinaryPowerBSSR(
    p = (th * q1 + q2) / (1 + th), Delta.A = 0.2, Delta.T = 0.2, N1 = ss$N1, N2 = ss$N2,
    n.interim = c(th * n12, n12), r = th, alpha = 0.05, tar.power = 0.8, Test = 'Chisq',
    alternative = 'two.sided', ss.method = 'null.variance', rounding = 'friede-kieser',
    N.max = 2 * ss$N
  )
  fk$n[i] <- ss$N
  fk$EN[i] <- res$E.N
  fk$power[i] <- res$power.BSSR
}
write.csv(fk, file.path(out.dir, 'friede-kieser-2004-table1.csv'), row.names = FALSE)
lab <- sprintf('FK2004 Table I: theta = %g, p1 = %g, n1 = %g', fk$theta, fk$p1, fk$n.pilot)
add('Friede and Kieser (2004)', paste0(lab, ': n'), fk$n.pub, fk$n, 0)
add('Friede and Kieser (2004)', paste0(lab, ': E[n]'), fk$EN.pub, fk$EN, 1)
add('Friede and Kieser (2004)', paste0(lab, ': power'), fk$pw.pub, fk$power, 3)

# Friede and Kieser (2004), Section 5 ---------------------------------------------------
# Depression trial: alternative 0.15, initial placebo rate 0.15, pilot of 120 patients.
# The article plots the level of the recalculation design against the required sample
# size n, with the upper bound 2n of the recalculated size; n is converted to the
# overall rate with formula (1)
plan <- BinarySampleSize(0.15 + 0.15, 0.15, 1, 0.05, 0.8, 'Chisq',
                         alternative = 'two.sided', method = 'null.variance')
add('Friede and Kieser (2004)', 'FK2004 Section 5: initial total sample size', 244,
    plan$N, 0)
rec <- BinaryBSSR(60, 60, round(0.467 * 120), 0.15, 1, 0.05, 0.8, 'Chisq',
                  alternative = 'two.sided', ss.method = 'null.variance',
                  rounding = 'friede-kieser')
add('Friede and Kieser (2004)', 'FK2004 Section 5: recalculated total sample size', 348,
    rec$N.final, 0,
    'Interim rate 56 / 120 = 0.4667, which the article reports as 0.467')
p.grid <- seq(0.010, 0.500, by = 0.001)
fixed.level <- max(BinaryPower(p.grid, p.grid, 122, 122, 0.05, 'Chisq',
                               alternative = 'two.sided')$Power)
add('Friede and Kieser (2004)', 'FK2004 Section 5: maximum level, fixed design', 0.0521,
    fixed.level, 4)
K <- (qnorm(0.975) + qnorm(0.8))^2
n.grid <- seq(20, 348, by = 2)
level_ipd <- function(level, level.ss) {
  vapply(n.grid, function(n) {
    v <- n * 0.15^2 / (4 * K)
    p <- (1 - sqrt(1 - 4 * v)) / 2
    BinaryPowerBSSR(p = p, Delta.A = 0.15, Delta.T = 0, N1 = 122, N2 = 122,
                    n.interim = c(60, 60), r = 1, alpha = level, tar.power = 0.8,
                    Test = 'Chisq', alternative = 'two.sided',
                    ss.method = 'null.variance', ss.alpha = level.ss,
                    rounding = 'friede-kieser', N.max = max(2 * n, 120))$power.BSSR
  }, numeric(1))
}
ipd <- level_ipd(0.05, 0.05)
# Level of the adjusted critical value 3.97 of the chi-squared statistic. The published
# maximum of 0.049 is reproduced when the recalculation also uses this level, as in
# Kieser (2020, Example 21.1); applying it to the final analysis only is reported as INFO
level.397 <- pchisq(3.97, 1, lower.tail = FALSE)
ipd.adj <- level_ipd(level.397, level.397)
ipd.adj.test <- level_ipd(level.397, 0.05)
write.csv(data.frame(n = n.grid, level = ipd, level.adjusted = ipd.adj,
                     level.adjusted.test.only = ipd.adj.test),
          file.path(out.dir, 'friede-kieser-2004-section5.csv'), row.names = FALSE)
add('Friede and Kieser (2004)', 'FK2004 Section 5: maximum level, internal pilot design',
    0.0518, max(ipd), 4)
add('Friede and Kieser (2004)', 'FK2004 Section 5: level of the critical value 3.97',
    0.046, level.397, 3)
add('Friede and Kieser (2004)',
    'FK2004 Section 5: maximum level with the critical value 3.97', 0.049, max(ipd.adj), 3,
    'The recalculation also uses the adjusted level')
add('Friede and Kieser (2004)',
    'FK2004 Section 5: maximum level with the critical value 3.97 in the final analysis only',
    0.049, max(ipd.adj.test), NA,
    'INFO: the recalculation uses the nominal level 0.05')

# Kieser (2020), Example 21.1 -------------------------------------------------------------
# BACLOREA trial: placebo rate 0.42 against 0.27, so group 1 (baclofen) has the lower
# rate and the alternative is 'less'. One-sided level 0.025, power 0.8, sample size from
# formula (21.3), pilot of 158 patients, unrestricted re-estimation
k.plan <- BinarySampleSize(0.27, 0.42, 1, 0.025, 0.8, 'Chisq', alternative = 'less',
                           method = 'standard')
add('Kieser (2020)', 'Kieser 2020: initial total sample size, Delta = -0.15', 314,
    k.plan$N, 0)
k.grid <- seq(0.01, 0.99, by = 0.005)
k.args <- list(Delta.A = -0.15, N1 = 157, N2 = 157, n.interim = c(79, 79), r = 1,
               alpha = 0.025, tar.power = 0.8, Test = 'Chisq', alternative = 'less',
               ss.method = 'standard', theta = k.grid)
tie <- do.call(BinaryTypeIErrorBSSR, c(k.args, list(maximize = 'grid')))
write.csv(tie, file.path(out.dir, 'kieser-2020-level-delta015.csv'), row.names = FALSE)
add('Kieser (2020)', 'Kieser 2020: fixed design, Delta = 0.15, maximum level', 0.0256,
    max(tie$TIE.TRAD), 4)
k.top.fixed <- tie$theta[abs(tie$TIE.TRAD - max(tie$TIE.TRAD)) < 1e-12]
add('Kieser (2020)', 'Kieser 2020: fixed design, Delta = 0.15, rate of the maximum', 0.37,
    min(k.top.fixed), 3,
    'The level is symmetric about 0.5, so the maximum is attained at 0.37 and 0.63')
add('Kieser (2020)', 'Kieser 2020: IPS design, Delta = 0.15, maximum level', 0.0267,
    max(tie$TIE.BSSR), 4)
k.top <- tie$theta[abs(tie$TIE.BSSR - max(tie$TIE.BSSR)) < 1e-12]
add('Kieser (2020)', 'Kieser 2020: IPS design, Delta = 0.15, rate of the maximum', 0.075,
    min(k.top), 3,
    'The level is symmetric about 0.5, so the maximum is attained at 0.075 and 0.925')
# The book reports the mean and the minimum over the overall rates that Figure 21.2
# shows, those for which the assumed rates of both groups lie in the unit interval:
# [0.075, 0.925] for Delta = 0.15 and [0.15, 0.85] for Delta = 0.30
in.range <- k.grid >= 0.075 - 1e-9 & k.grid <= 0.925 + 1e-9
range.note <- 'Over the overall rates shown in Figure 21.2 a, [0.075, 0.925]'
add('Kieser (2020)', 'Kieser 2020: fixed design, Delta = 0.15, mean level', 0.0248,
    mean(tie$TIE.TRAD[in.range]), 4, range.note)
add('Kieser (2020)', 'Kieser 2020: IPS design, Delta = 0.15, mean level', 0.0250,
    mean(tie$TIE.BSSR[in.range]), 4, range.note)
add('Kieser (2020)', 'Kieser 2020: fixed design, Delta = 0.15, minimum level', 0.0236,
    min(tie$TIE.TRAD[in.range]), 4, range.note)
add('Kieser (2020)', 'Kieser 2020: IPS design, Delta = 0.15, minimum level', 0.0235,
    min(tie$TIE.BSSR[in.range]), 4, range.note)
for (adj in c('test', 'both')) {
  a <- do.call(BinaryAlphaAdjBSSR, c(k.args, list(adjust = adj, step = 1e-5)))
  write.csv(a, file.path(out.dir, sprintf('kieser-2020-adjusted-%s.csv', adj)),
            row.names = FALSE)
  if (adj == 'test') {
    add('Kieser (2020)', 'Kieser 2020: adjusted level, fixed design', 0.0243,
        a$alpha.adj[a$Design == 'TRAD'], 4)
  }
  # The book does not state whether the recalculation uses the adjusted level. The value
  # it reports is compared with adjust = 'both', and the other variant is reported as INFO
  add('Kieser (2020)', sprintf('Kieser 2020: adjusted level, IPS design (adjust = %s)', adj),
      0.0247, a$alpha.adj[a$Design == 'BSSR'], if (adj == 'both') 4 else NA,
      if (adj == 'both') {
        'The recalculation also uses the adjusted level'
      } else {
        'INFO: the adjusted level applied to the final analysis only'
      })
}
# Same design with Delta = 0.30, n0 = 2 x 39 and a pilot of 40 patients
k3 <- do.call(BinaryTypeIErrorBSSR,
              modifyList(k.args, list(Delta.A = -0.30, N1 = 39, N2 = 39,
                                      n.interim = c(20, 20), maximize = 'grid')))
write.csv(k3, file.path(out.dir, 'kieser-2020-level-delta030.csv'), row.names = FALSE)
mid <- k.grid >= 0.15 - 1e-9 & k.grid <= 0.85 + 1e-9
range.note <- 'Over the overall rates shown in Figure 21.2 b, [0.15, 0.85]'
add('Kieser (2020)', 'Kieser 2020: fixed design, Delta = 0.30, maximum level', 0.0268,
    max(k3$TIE.TRAD[mid]), 4, range.note)
add('Kieser (2020)', 'Kieser 2020: fixed design, Delta = 0.30, minimum level', 0.0253,
    min(k3$TIE.TRAD[mid]), 4, range.note)
add('Kieser (2020)', 'Kieser 2020: fixed design, Delta = 0.30, mean level', 0.0260,
    mean(k3$TIE.TRAD[mid]), 4, range.note)
add('Kieser (2020)', 'Kieser 2020: IPS design, Delta = 0.30, maximum level', 0.0273,
    max(k3$TIE.BSSR[mid]), 4, range.note)
add('Kieser (2020)', 'Kieser 2020: IPS design, Delta = 0.30, minimum level', 0.0259,
    min(k3$TIE.BSSR[mid]), 4, range.note)
add('Kieser (2020)', 'Kieser 2020: IPS design, Delta = 0.30, mean level', 0.0268,
    mean(k3$TIE.BSSR[mid]), 4, range.note)
add('Kieser (2020)', 'Kieser 2020: fixed design, Delta = 0.30, maximum over [0.01, 0.99]',
    NA, max(k3$TIE.TRAD), NA, 'INFO')
# Figure 21.1: actual levels at the assumed overall rate for the alternatives 0.15 to 0.30
# (pE > pC), the assumed overall rates 0.30 to 0.50 and pilots of 25, 50 and 75 per cent
# of the fixed sample size, with ceiling(fraction nC) patients in group C and r times as
# many in group E, as in inst/reproduce/reproduce-figures.R. The book reports the mean,
# minimum and maximum over these designs. The book does not state how the interim
# outcomes with a negative recovered rate of group C are treated: bbssr truncates the
# Bernoulli variance of that group at zero, and the levels when the trial stops after the
# pilot for such outcomes are reported as INFO
k1.stop <- function(p, Delta, r, N1, N2, n.interim) {
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
k1 <- expand.grid(pA = round(seq(0.30, 0.50, by = 0.01), 2),
                  Delta = c(0.15, 0.20, 0.25, 0.30), fraction = c(0.25, 0.5, 0.75),
                  r = c(1, 3))
k1[c('fixed', 'ips', 'ips.stop')] <- NA_real_
for (i in seq_len(nrow(k1))) {
  ri <- k1$r[i]
  pE <- k1$pA[i] + k1$Delta[i] / (1 + ri)
  pC <- k1$pA[i] - ri * k1$Delta[i] / (1 + ri)
  ss <- BinarySampleSize(pE, pC, ri, 0.025, 0.8, 'Chisq', method = 'standard')
  n1C <- ceiling(k1$fraction[i] * ss$N2)
  res <- BinaryPowerBSSR(p = k1$pA[i], Delta.A = k1$Delta[i], Delta.T = 0, N1 = ss$N1,
                         N2 = ss$N2, n.interim = c(ri * n1C, n1C), r = ri, alpha = 0.025,
                         tar.power = 0.8, Test = 'Chisq', ss.method = 'standard')
  k1[i, c('fixed', 'ips', 'ips.stop')] <- c(
    res$power.TRAD, res$power.BSSR,
    k1.stop(k1$pA[i], k1$Delta[i], ri, ss$N1, ss$N2, c(ri * n1C, n1C))
  )
}
write.csv(k1, file.path(out.dir, 'kieser-2020-figure21-1.csv'), row.names = FALSE)
k1.add <- function(label, published, x, digits = 4, note = '') {
  add('Kieser (2020)',
      sprintf('Kieser 2020 Figure 21.1: %s, %s level', label,
              c('mean', 'minimum', 'maximum')),
      published, c(mean(x), min(x), max(x)), digits, note)
}
# The fixed design does not depend on the pilot, so it is summarized over one fraction
k1.add('r = 1, fixed design', c(0.0257, 0.0233, 0.0308),
       k1$fixed[k1$r == 1 & k1$fraction == 0.5])
k1.add('r = 1, IPS design', c(0.0257, 0.0235, 0.0289), k1$ips[k1$r == 1])
k1.pub <- list(c(0.0256, 0.0238, 0.0279), c(0.0257, 0.0236, 0.0286),
               c(0.0257, 0.0235, 0.0289))
for (j in 1:3) {
  fr <- c(0.25, 0.5, 0.75)[j]
  k1.add(sprintf('r = 1, IPS design, pilot %g%%', 100 * fr), k1.pub[[j]],
         k1$ips[k1$r == 1 & k1$fraction == fr])
}
k1.add('r = 3, fixed design', c(0.0243, 0.0214, 0.0279),
       k1$fixed[k1$r == 3 & k1$fraction == 0.5])
k1.add('r = 3, IPS design', c(0.0241, 0.0174, 0.0274), k1$ips[k1$r == 3])
for (ri in c(1, 3)) {
  k1.add(sprintf('r = %d, IPS design stopped after the pilot', ri),
         if (ri == 1) c(0.0257, 0.0235, 0.0289) else c(0.0241, 0.0174, 0.0274),
         k1$ips.stop[k1$r == ri], NA,
         'INFO: the trial stops after the pilot when a recovered rate is negative')
}
k1.reason <- paste(
  'The book does not state how the interim outcomes with a negative recovered rate of',
  'group C are treated. bbssr truncates the Bernoulli variance of that group at zero and',
  're-estimates the sample size. If the trial stops after the pilot for such outcomes,',
  'the minimum of the book is reproduced (INFO items of the same figure and the vignette',
  'Validation by Reproducing Published Figures).')
explain('Kieser 2020 Figure 21.1: r = 3, IPS design, minimum level', 0.0035, k1.reason)
explain('Kieser 2020 Figure 21.1: r = 3, IPS design, mean level', 1e-4, k1.reason)
# Distribution of the recalculated total sample size (Figure 21.4). The book reports
# totals such as 315, so the total is rounded up rather than each group
kN <- BinaryPowerBSSR(p = c(0.25, 0.35, 0.45), Delta.A = -0.15, Delta.T = -0.15,
                      N1 = 157, N2 = 157, n.interim = c(79, 79), r = 1, alpha = 0.025,
                      tar.power = 0.8, Test = 'Chisq', alternative = 'less',
                      ss.method = 'standard', rounding = 'total')
kS <- summary(kN)
write.csv(kS, file.path(out.dir, 'kieser-2020-sample-size.csv'), row.names = FALSE)
add('Kieser (2020)', sprintf('Kieser 2020: median total sample size, p = %g', kS$p),
    c(258, 315, 343), kS$N.q50, 0)
add('Kieser (2020)', 'Kieser 2020: lower quartile of the total sample size, p = 0.35', 303,
    kS$N.q25[2], 0)
add('Kieser (2020)', 'Kieser 2020: upper quartile of the total sample size, p = 0.35', 325,
    kS$N.q75[2], 0)
add('Kieser (2020)', sprintf('Kieser 2020: interquartile range, p = %g', kS$p[c(1, 3)]),
    c(31, 7), kS$N.q75[c(1, 3)] - kS$N.q25[c(1, 3)], 0)

# Farrington and Manning (1990), Tables I and II and the first example -------------------
# One-sided level 0.05 and power 0.9 for the risk difference s0 under the null hypothesis,
# with sample sizes from formula (4) rounded to the nearest integer. The article tests
# p1 - p2 <= s0 and writes theta = N2 / N1, so the margin is -s0 and r = 1 / theta. The
# true powers are computed at the published sample sizes. The published powers agree with
# the powers without the outcome of no responders in either group, which the test rejects
# in the settings with s0 < 0, truncated to two decimals (in per cent), except for the
# entry p1 = 0.5, p2 = 0.1, s0 = 0.2, theta = 3/2, which lies 0.01 percentage points below
# the recomputed power. The recomputed power leaves this outcome out and is compared after
# truncation
fm.src <- 'Farrington and Manning (1990)'
fm <- data.frame(
  p1 = rep(c(0.1, 0.2, 0.5, 0.05, 0.1, 0.25, 0.01, 0.02, 0.05), each = 3),
  p2 = rep(c(0.1, 0.1, 0.1, 0.05, 0.05, 0.05, 0.01, 0.01, 0.01), each = 3),
  s0 = rep(c(-0.2, -0.1, 0.2, -0.2, -0.05, 0.1, -0.2, -0.01, 0.02), each = 3),
  theta = rep(c(2 / 3, 1, 3 / 2), times = 9),
  N1.pub = c(63, 46, 34, 75, 57, 45, 94, 78, 67, 45, 32, 23, 168, 128, 101, 235, 197, 173,
             26, 18, 13, 914, 695, 544, 1357, 1151, 1018),
  N2.pub = c(42, 46, 52, 50, 57, 68, 62, 78, 101, 30, 32, 35, 112, 128, 151, 156, 197, 260,
             18, 18, 19, 609, 695, 816, 905, 1151, 1527),
  pw.pub = c(91.56, 91.17, 90.98, 91.25, 91.34, 91.20, 90.14, 90.60, 90.23, 88.40, 88.84,
             90.70, 90.84, 91.36, 92.68, 90.77, 90.59, 90.65, 22.68, 16.34, 12.07, 90.75,
             91.97, 92.63, 90.91, 90.88, 90.68) / 100
)
fm[c('N1', 'N2', 'power', 'power.article')] <- NA_real_
for (i in seq_len(nrow(fm))) {
  m <- -fm$s0[i]
  ss <- BinarySampleSize(fm$p1[i], fm$p2[i], 1 / fm$theta[i], 0.05, 0.9,
                         'Farrington-Manning', method = 'standard', rounding = 'nearest',
                         margin = m)
  fm$N1[i] <- ss$N1
  fm$N2[i] <- ss$N2
  n1 <- fm$N1.pub[i]
  n2 <- fm$N2.pub[i]
  fm$power[i] <- BinaryPower(fm$p1[i], fm$p2[i], n1, n2, 0.05, 'Farrington-Manning',
                             margin = m)$Power
  # The first cell of the rejection region is the outcome (0, 0)
  rej00 <- as.vector(BinaryRR(n1, n2, 0.05, 'Farrington-Manning', margin = m))[1]
  fm$power.article[i] <- fm$power[i] -
    rej00 * dbinom(0, n1, fm$p1[i]) * dbinom(0, n2, fm$p2[i])
}
write.csv(fm, file.path(out.dir, 'farrington-manning-1990-table1.csv'), row.names = FALSE)
fm.lab <- sprintf('FM1990 Table I: p1 = %g, p2 = %g, s0 = %g, theta = %s', fm$p1, fm$p2,
                  fm$s0, rep(c('2/3', '1', '3/2'), times = 9))
add(fm.src, paste0(fm.lab, ': N1'), fm$N1.pub, fm$N1, 0)
add(fm.src, paste0(fm.lab, ': N2'), fm$N2.pub, fm$N2, 0)
fm.note <- ifelse(fm$power - fm$power.article > 5e-5,
                  sprintf(paste('Without the outcome (0, 0); bbssr, which rejects it,',
                                'gives %.4f'), fm$power), '')
# The entry that is not reproduced stays a FAIL, since no reason has been found
fm.note[fm$p1 == 0.5 & fm$p2 == 0.1 & abs(fm$theta - 1.5) < 1e-12] <- paste(
  'Not reproduced, cause not identified. An independent implementation with numpy and',
  'scipy gives the same 0.90241. Neighbouring sample sizes (66 to 68 and 100 to 102),',
  'the critical value 1.645, exchanged groups and single-precision arithmetic do not',
  'give the published value.')
add(fm.src, paste0(fm.lab, ': true power'), fm$pw.pub, fm$power.article, 4, fm.note,
    rule = 'truncate')
explain(paste0(fm.lab[c(5, 20)], ': true power'), 2e-5, paste(
  'The recomputed power lies less than 0.002 percentage points below the published',
  'value, so truncation to two decimals gives one unit less; rounding gives the',
  'published value.'))
# First example: p1 = 0.4, p2 = 0.05 and s0 = 0.2 with equal groups, level 0.05 and
# power 0.8
fm.ex <- bbssr:::fm_restricted(0.4, 0.05, 1, 0.2)
add(fm.src, 'FM1990 example: restricted estimate of p1', 0.2935, fm.ex$p1, 4)
add(fm.src, 'FM1990 example: restricted estimate of p2', 0.0935, fm.ex$p2, 4)
ss <- BinarySampleSize(0.4, 0.05, 1, 0.05, 0.8, 'Farrington-Manning', method = 'standard',
                       rounding = 'nearest', margin = -0.2)
add(fm.src, c('FM1990 example: N1', 'FM1990 example: N2'), 80, c(ss$N1, ss$N2), 0)
add(fm.src, 'FM1990 example: true power at N1 = N2 = 80', 0.813,
    BinaryPower(0.4, 0.05, 80, 80, 0.05, 'Farrington-Manning', margin = -0.2)$Power, 3)
# Table II, Method 3: the two settings with equal groups that Table I does not contain
fm2 <- data.frame(p = c(0.05, 0.01), s0 = c(-0.1, -0.02), N.pub = c(103, 558))
fm2[c('N1', 'N2')] <- NA_real_
for (i in seq_len(nrow(fm2))) {
  ss <- BinarySampleSize(fm2$p[i], fm2$p[i], 1, 0.05, 0.9, 'Farrington-Manning',
                         method = 'standard', rounding = 'nearest', margin = -fm2$s0[i])
  fm2$N1[i] <- ss$N1
  fm2$N2[i] <- ss$N2
}
fm2.lab <- sprintf('FM1990 Table II: p1 = p2 = %g, s0 = %g, Method 3', fm2$p, fm2$s0)
add(fm.src, paste0(fm2.lab, ': N1'), fm2$N.pub, fm2$N1, 0)
add(fm.src, paste0(fm2.lab, ': N2'), fm2$N.pub, fm2$N2, 0)
# Table II, Method 1 for the difference (p. 1450): the null variance is also evaluated at
# the true values p1 and p2, which is the formula of method = 'alternative.variance', with
# equal groups and each group rounded to the nearest whole number. Method 2 (fixed
# marginal totals) is not implemented: the article does not recommend it, and it fails
# for three of the nine settings of the table
fm3 <- data.frame(p1 = c(0.1, 0.2, 0.5, 0.05, 0.1, 0.25, 0.01, 0.02, 0.05),
                  p2 = c(0.1, 0.1, 0.1, 0.05, 0.05, 0.05, 0.01, 0.01, 0.01),
                  s0 = c(-0.2, -0.1, 0.2, -0.1, -0.05, 0.1, -0.02, -0.01, 0.02),
                  N.pub = c(39, 54, 73, 81, 118, 201, 424, 632, 1229))
fm3[c('N1', 'N2')] <- NA_real_
for (i in seq_len(nrow(fm3))) {
  ss <- BinarySampleSize(fm3$p1[i], fm3$p2[i], 1, 0.05, 0.9, 'Blackwelder',
                         method = 'alternative.variance', rounding = 'nearest',
                         margin = -fm3$s0[i])
  fm3$N1[i] <- ss$N1
  fm3$N2[i] <- ss$N2
}
fm3.lab <- sprintf('FM1990 Table II: p1 = %g, p2 = %g, s0 = %g, Method 1', fm3$p1, fm3$p2,
                   fm3$s0)
add(fm.src, paste0(fm3.lab, ': N1'), fm3$N.pub, fm3$N1, 0)
add(fm.src, paste0(fm3.lab, ': N2'), fm3$N.pub, fm3$N2, 0)

# Blackwelder (1982), Table 3 and the examples ------------------------------------------
# The statistics use the unpooled standard error. In the examples with 30 patients per
# group (pp. 347 and 350), the rows and columns of the p-value matrix are the numbers of
# responders plus one. Group 1 is the standard therapy for the conventional null
# hypothesis H0 and the experimental therapy for the null hypothesis H0' that the standard
# therapy is better by at least delta
bw.src <- 'Blackwelder (1982)'
p.sup <- attr(BinaryRR(30, 30, 0.05, 'Blackwelder'), 'p.value')
p.ni <- attr(BinaryRR(30, 30, 0.05, 'Blackwelder', margin = 0.2), 'p.value')
add(bw.src, 'BW1982 p. 347: statistic for H0, 18 / 30 against 13 / 30', 1.31,
    qnorm(p.sup[19, 14], lower.tail = FALSE), 2)
add(bw.src, 'BW1982 p. 347: one-sided p-value for H0, 18 / 30 against 13 / 30', 0.095,
    p.sup[19, 14], 3, 'The article gives about 0.095')
add(bw.src, 'BW1982 p. 350: statistic for H0, 21 / 30 against 18 / 30', 0.816,
    qnorm(p.sup[22, 19], lower.tail = FALSE), 3)
add(bw.src, 'BW1982 p. 350: one-sided p-value for H0, 21 / 30 against 18 / 30', 0.21,
    p.sup[22, 19], 2)
add(bw.src, "BW1982 p. 350: statistic for H0' with delta = 0.2, 18 / 30 against 21 / 30",
    -0.816, -qnorm(p.ni[19, 22], lower.tail = FALSE), 3,
    paste('The article subtracts the experimental rate from the standard rate, the',
          'reverse of bbssr'))
add(bw.src, paste("BW1982 p. 350: one-sided p-value for H0' with delta = 0.2, 18 / 30",
                  'against 21 / 30'), 0.21, p.ni[19, 22], 2)
# Table 3: one-sided level 0.05, power 0.9 and equal groups, each rounded up. The true
# difference is delta under H0 and ps - pe under H0'
bw.sup <- data.frame(ps = c(0.9, 0.9, 0.6, 0.6, 0.4, 0.4),
                     pe = c(0.8, 0.7, 0.5, 0.4, 0.3, 0.2),
                     N.pub = c(430, 130, 840, 206, 772, 172))
bw.sup$N <- vapply(seq_len(nrow(bw.sup)), function(i) {
  BinarySampleSize(bw.sup$ps[i], bw.sup$pe[i], 1, 0.05, 0.9, 'Blackwelder',
                   method = 'alternative.variance')$N
}, numeric(1))
bw.ni <- data.frame(ps = c(0.9, 0.9, 0.6, 0.6, 0.4, 0.4, 0.9, 0.9, 0.6, 0.6),
                    pe = c(0.9, 0.9, 0.6, 0.6, 0.4, 0.4, 0.85, 0.8, 0.55, 0.5),
                    delta = c(0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2),
                    N.pub = c(310, 78, 824, 206, 824, 206, 1492, 430, 3342, 840))
bw.ni$N <- vapply(seq_len(nrow(bw.ni)), function(i) {
  BinarySampleSize(bw.ni$pe[i], bw.ni$ps[i], 1, 0.05, 0.9, 'Blackwelder',
                   method = 'alternative.variance', margin = bw.ni$delta[i])$N
}, numeric(1))
write.csv(rbind(data.frame(hypothesis = 'H0', bw.sup[c('ps', 'pe')], delta = 0,
                           bw.sup[c('N.pub', 'N')]),
                data.frame(hypothesis = "H0'", bw.ni)),
          file.path(out.dir, 'blackwelder-1982-table3.csv'), row.names = FALSE)
add(bw.src, sprintf('BW1982 Table 3, H0: ps = %g, pe = %g: total sample size', bw.sup$ps,
                    bw.sup$pe), bw.sup$N.pub, bw.sup$N, 0)
bw.lab <- sprintf("BW1982 Table 3, H0': ps = %g, pe = %g, delta = %g: total sample size",
                  bw.ni$ps, bw.ni$pe, bw.ni$delta)
add(bw.src, bw.lab, bw.ni$N.pub, bw.ni$N, 0)
add(bw.src, paste("BW1982 Table 3, H0': ps = 0.6, pe = 0.55, delta = 0.1: formula of the",
                  'article with the quantiles 1.645 and 1.282'), 3342,
    2 * ceiling((1.645 + 1.282)^2 * (0.6 * 0.4 + 0.55 * 0.45) / (0.6 - 0.55 - 0.1)^2), 0,
    'Evaluated in this script, not by bbssr')
explain(bw.lab[9], 2, paste(
  'The article uses the rounded quantiles 1.645 and 1.282, with which its formula gives',
  '3342 (next item); the exact quantiles give 3340.'))

# Friede, Mitchell and Mueller-Velten (2007), Tables 2 and 3 and Sections 5 and 6 --------
# Margin 0.1, one-sided level 0.025, target power 0.8 and assumed difference 0, with the
# overall response rate as the nuisance parameter. Group 1 is the experimental group, and
# the article writes q = n2 / n1, so r = 1 / q. The fixed sample sizes round the two groups
# up separately (rounding = 'friede-kieser'). The article does not state how the
# re-estimated sample sizes are rounded; here each group is rounded to the nearest integer
# (rounding = 'nearest')
fmm.src <- 'Friede et al. (2007)'
fmm_n <- function(p, q, Test, method) {
  BinarySampleSize(p, p, 1 / q, 0.025, 0.8, Test, method = method,
                   rounding = 'friede-kieser', margin = 0.1)
}
fmm_bssr <- function(p, N1, N2, n.interim, r, Test, ss.method) {
  BinaryPowerBSSR(p = p, Delta.A = 0, Delta.T = 0, N1 = N1, N2 = N2, n.interim = n.interim,
                  r = r, alpha = 0.025, tar.power = 0.8, Test = Test,
                  ss.method = ss.method, rounding = 'nearest', margin = 0.1)
}
# Section 5: sample sizes of formulae (1) (Blackwelder) and (2) (Farrington-Manning) quoted
# in the text
s5 <- data.frame(p = c(0.7, 0.7, 0.9, 0.9, 0.9, 0.9), q = c(1, 1 / 3, 1, 1 / 3, 1, 1 / 3),
                 Test = rep(c('Blackwelder', 'Farrington-Manning'), times = c(4, 2)),
                 N.pub = c(660, 880, 284, 378, 310, 306), stringsAsFactors = FALSE)
s5$N <- vapply(seq_len(nrow(s5)), function(i) {
  fmm_n(s5$p[i], s5$q[i], s5$Test[i],
        if (s5$Test[i] == 'Blackwelder') 'alternative.variance' else 'standard')$N
}, numeric(1))
add(fmm.src, sprintf('FMM2007 Section 5: %s test, p = %g, q = %s: total sample size',
                     s5$Test, s5$p, rep(c('1', '1/3'), times = 3)), s5$N.pub, s5$N, 0)
# Table 2: Farrington-Manning test, true overall rate equal to the assumed rate and an
# internal pilot study of half of each group of the fixed design, rounded to the nearest
# integer
t2 <- data.frame(pa = rep(c(0.3, 0.4, 0.5, 0.6, 0.7), each = 3),
                 q = rep(c(1 / 3, 1 / 2, 1), times = 5),
                 pw.fixed.pub = c(80.1, 80.2, 80.1, 80.3, 80.1, 80.4, 80.1, 79.9, 79.5,
                                  80.2, 80.0, 80.4, 80.2, 80.1, 80.1) / 100,
                 pw.bssr.pub = c(rep(80.0, 8), 79.5, rep(80.0, 6)) / 100,
                 N.pub = c(928, 770, 658, 1023, 857, 750, 1035, 876, 780, 966, 825, 750,
                           815, 707, 658),
                 q05.pub = c(875, 717, 602, 997, 831, 720, 1018, 863, 772, 920, 786, 720,
                             736, 641, 602),
                 EN.pub = c(925.2, 767.3, 655.8, 1019.5, 854.5, 747.0, 1031.8, 872.9, 776.9,
                            962.2, 822.5, 747.0, 811.7, 703.9, 655.8),
                 q95.pub = c(971, 812, 704, 1035, 872, 768, 1039, 877, 778, 999, 852, 768,
                             881, 761, 704))
t2[c('N', 'pw.fixed', 'pw.bssr', 'q05', 'EN', 'q95')] <- NA_real_
for (i in seq_len(nrow(t2))) {
  ss <- fmm_n(t2$pa[i], t2$q[i], 'Farrington-Manning', 'standard')
  res <- fmm_bssr(t2$pa[i], ss$N1, ss$N2, floor(c(ss$N1, ss$N2) / 2 + 0.5), 1 / t2$q[i],
                  'Farrington-Manning', 'standard')
  sm <- summary(res, probs = c(0.05, 0.95))
  qs <- grep('^N\\.q', names(sm))
  t2$N[i] <- ss$N
  t2$pw.fixed[i] <- res$power.TRAD
  t2$pw.bssr[i] <- res$power.BSSR
  t2$q05[i] <- sm[1, qs[1]]
  t2$EN[i] <- res$E.N
  t2$q95[i] <- sm[1, qs[2]]
}
write.csv(t2, file.path(out.dir, 'friede-2007-table2.csv'), row.names = FALSE)
t2.lab <- sprintf('FMM2007 Table 2: pa = %g, q = %s', t2$pa, rep(c('1/3', '1/2', '1'), 5))
t2.it <- function(x) paste0(t2.lab, ': ', x)
add(fmm.src, t2.it('power, fixed design'), t2$pw.fixed.pub, t2$pw.fixed, 3)
add(fmm.src, t2.it('power, re-estimation'), t2$pw.bssr.pub, t2$pw.bssr, 3)
add(fmm.src, t2.it('total sample size, fixed design'), t2$N.pub, t2$N, 0)
add(fmm.src, t2.it('5% quantile of the total sample size'), t2$q05.pub, t2$q05, 0)
add(fmm.src, t2.it('mean total sample size'), t2$EN.pub, t2$EN, 1)
add(fmm.src, t2.it('95% quantile of the total sample size'), t2$q95.pub, t2$q95, 0)
fmm.round <- paste(
  'The article does not state how the re-estimated sample sizes are rounded to whole',
  'patients. The rounding rule shifts the distribution of the total sample size by one or',
  'two patients and the power by fractions of a percentage point; here each group is',
  'rounded to the nearest integer.')
explain(t2.it('power, re-estimation'), 1e-3, fmm.round)
explain(t2.it('mean total sample size'), 0.7, fmm.round)
explain(c(t2.it('5% quantile of the total sample size'),
          t2.it('95% quantile of the total sample size')), 2, fmm.round)
# Table 3: Farrington-Manning test, assumed overall rate 0.7 and true overall rate 0.42,
# equal groups of n0 / 2 patients and an internal pilot study of n0 / 5 patients
t3 <- data.frame(n0 = seq(400, 1000, by = 100),
                 pw.fixed.pub = c(52.34, 62.46, 70.46, 76.67, 81.65, 86.46, 89.44) / 100,
                 pw.bssr.pub = c(79.45, 79.55, 79.62, 79.66, 79.68, 79.72, 79.75) / 100,
                 EN.pub = c(750.7, 752.6, 753.6, 754.4, 755.2, 755.7, 756.1),
                 SD.pub = c(29.5, 26.1, 23.4, 21.6, 20.1, 18.8, 17.6))
t3[c('pw.fixed', 'pw.bssr', 'EN', 'SD', 'N.max')] <- NA_real_
for (i in seq_len(nrow(t3))) {
  res <- fmm_bssr(0.42, t3$n0[i] / 2, t3$n0[i] / 2, rep(t3$n0[i] / 10, 2), 1,
                  'Farrington-Manning', 'standard')
  map <- attr(res, 'reestimation')
  t3$pw.fixed[i] <- res$power.TRAD
  t3$pw.bssr[i] <- res$power.BSSR
  t3$EN[i] <- res$E.N
  t3$SD[i] <- summary(res)$SD.N
  t3$N.max[i] <- max(map$N1 + map$N2)
}
write.csv(t3, file.path(out.dir, 'friede-2007-table3.csv'), row.names = FALSE)
t3.it <- function(x) paste0(sprintf('FMM2007 Table 3: n0 = %g', t3$n0), ': ', x)
add(fmm.src, t3.it('power, fixed design'), t3$pw.fixed.pub, t3$pw.fixed, 4)
add(fmm.src, t3.it('power, re-estimation'), t3$pw.bssr.pub, t3$pw.bssr, 4)
add(fmm.src, t3.it('mean total sample size'), t3$EN.pub, t3$EN, 1)
add(fmm.src, t3.it('standard deviation of the total sample size'), t3$SD.pub, t3$SD, 1)
add(fmm.src, t3.it('largest total sample size'), 780, t3$N.max, 0,
    'The caption gives 780 for all lines')
fmm.bound <- paste(
  'The article computes the power only to within 0.0001, as the mean of an upper and a',
  'lower bound (Section 4).')
explain(t3.it('power, fixed design'), 1e-4, fmm.bound)
explain(t3.it('power, re-estimation'), 2e-4, paste(fmm.bound, fmm.round))
explain(t3.it('mean total sample size'), 0.4, fmm.round)
explain(t3.it('standard deviation of the total sample size'), 0.2, fmm.round)
# Section 6: the example trial with 330 patients per group, true overall rate 0.42 and an
# internal pilot study of 20 per cent (66 patients per group). The Blackwelder test
# re-estimates with formula (1) and the Farrington-Manning test with formula (2)
s6.pub <- list('Blackwelder' = c(0.745, 0.797, 759),
               'Farrington-Manning' = c(0.749, 0.797, 754))
for (tst in names(s6.pub)) {
  res <- fmm_bssr(0.42, 330, 330, c(66, 66), 1, tst,
                  if (tst == 'Blackwelder') 'alternative.variance' else 'standard')
  s6.lab <- sprintf('FMM2007 Section 6: %s test: %s', tst,
                    c('power, fixed design', 'power, re-estimation',
                      'expected total sample size'))
  add(fmm.src, s6.lab, s6.pub[[tst]], c(res$power.TRAD, res$power.BSSR, res$E.N),
      c(3, 3, 0))
}

# Boschloo (1970), Sections 2, 3, 4 and 6 ------------------------------------------------
# Fisher's test is used at a raised conditional level gamma, the largest level at which
# the unconditional size does not exceed alpha. Rejecting when the Fisher p-value is at
# most gamma is the same as rejecting when the Boschloo p-value is at most alpha. The
# article writes m and n for the sizes of the first and the second sample and tests
# p1 <= p2 against p1 > p2, so N1 = m, N2 = n and the alternative is 'greater'. The
# two-sided test gives each of the two parts of the critical region the conditional level
# gamma / 2, which is the 'central' convention of tsmethod. The rows and columns of the
# matrices of bbssr are the numbers of successes plus one
bo.src <- 'Boschloo (1970)'
bo_cells <- function(RR) matrix(as.vector(RR), attr(RR, 'N1') + 1L, attr(RR, 'N2') + 1L)
# Section 4: power with m = n = 15 of Fisher's test (column I) and of Fisher's test at the
# raised level (column II). Column III, the randomized test, is not part of the package
bo <- data.frame(alpha = rep(c(0.01, 0.05), each = 4), p1 = rep(c(0.3, 0.6, 0.7, 0.8), 2),
                 p2 = rep(c(0.1, 0.1, 0.2, 0.2), 2),
                 fisher.pub = c(0.0755, 0.6087, 0.5268, 0.7647,
                                0.2558, 0.8451, 0.8066, 0.9391),
                 boschloo.pub = c(0.1493, 0.7394, 0.6796, 0.8723,
                                  0.3531, 0.9189, 0.8962, 0.9744))
bo[c('fisher', 'boschloo')] <- NA_real_
for (a in c(0.01, 0.05)) {
  i <- bo$alpha == a
  bo$fisher[i] <- BinaryPower(bo$p1[i], bo$p2[i], 15, 15, a, 'Fisher')$Power
  bo$boschloo[i] <- BinaryPower(bo$p1[i], bo$p2[i], 15, 15, a, 'Boschloo')$Power
}
write.csv(bo, file.path(out.dir, 'boschloo-1970-section4.csv'), row.names = FALSE)
bo.lab <- sprintf('Boschloo1970 Section 4: alpha = %g, p1 = %g, p2 = %g', bo$alpha, bo$p1,
                  bo$p2)
add(bo.src, paste0(bo.lab, ': power, Fisher test (column I)'), bo$fisher.pub, bo$fisher,
    4)
add(bo.src, paste0(bo.lab, ': power, raised level (column II)'), bo$boschloo.pub,
    bo$boschloo, 4)
explain(paste0(bo.lab[6], ': power, Fisher test (column I)'), 6e-5, paste(
  'The recomputed power lies less than 1e-7 above the rounding boundary 0.84515, so',
  'rounding to four decimals gives one unit more than the published value. The table',
  'rounds rather than truncates, as the entry 0.6087 shows (the power is 0.60866).'))
# Sections 2, 3 and 6: m = 15, n = 10 and the one-sided level 0.05. The example of Section
# 6 observes 5 successes out of 15 against 0 out of 10
bo.f <- BinaryRR(15, 10, 0.05, 'Fisher')
bo.b <- BinaryRR(15, 10, 0.05, 'Boschloo')
bo.f2 <- BinaryRR(15, 10, 0.05, 'Fisher', alternative = 'two.sided', tsmethod = 'central')
bo.b2 <- BinaryRR(15, 10, 0.05, 'Boschloo', alternative = 'two.sided',
                  tsmethod = 'central')
bo.pf <- attr(bo.f, 'p.value')
bo.pf2 <- attr(bo.f2, 'p.value')
bo.theta <- seq(0, 1, by = 1e-4)
bo.size <- max(vapply(bo.theta, function(t) {
  sum(outer(dbinom(0:15, 15, t), dbinom(0:10, 10, t)) * bo_cells(bo.f))
}, numeric(1)))
bo.ex <- 'Boschloo1970 Section 6: 5 / 15 against 0 / 10'
add(bo.src, paste0(bo.ex, ': Fisher p-value'), 0.0565, bo.pf[6, 1], 4)
add(bo.src, paste0(bo.ex, ': rejected by the Fisher test at 0.05 (1 = yes)'), 0,
    as.numeric(bo_cells(bo.f)[6, 1]), 0)
add(bo.src, paste0(bo.ex, ': rejected at the raised level, one-sided (1 = yes)'), 1,
    as.numeric(bo_cells(bo.b)[6, 1]), 0)
add(bo.src, paste0(bo.ex, ': rejected at the raised level, two-sided (1 = yes)'), 1,
    as.numeric(bo_cells(bo.b2)[6, 1]), 0, "tsmethod = 'central'")
add(bo.src, 'Boschloo1970 Section 2: size of the Fisher test, m = 15, n = 10, level 0.05',
    0.02, bo.size, 2, 'The article gives about .02; maximum over a grid of step 1e-4')
bo.fig <- 'Figure 1, m = 15, n = 10, one-sided level 0.05'
add(bo.src, paste('Boschloo1970 Section 6: outcomes added to the critical region of the',
                  'Fisher test'), 8, sum(bo_cells(bo.b) & !bo_cells(bo.f)), 0, bo.fig)
add(bo.src, paste('Boschloo1970 Section 6: outcomes removed from the critical region of',
                  'the Fisher test'), 0, sum(bo_cells(bo.f) & !bo_cells(bo.b)), 0, bo.fig)
# Raised levels quoted from the table of the article: 0.09 for the one-sided and 0.114 for
# the two-sided test at the level 0.05. The Fisher test at these conditional levels and
# the Boschloo test must have the same critical region
add(bo.src, paste('Boschloo1970 Sections 3 and 6: outcomes on which the Fisher test at',
                  'the raised level 0.09 and the Boschloo test differ, one-sided'),
    0, sum((bo.pf <= 0.09) != bo_cells(bo.b)), 0)
add(bo.src, paste('Boschloo1970 Section 6: outcomes on which the Fisher test at the',
                  'raised level 0.114 and the Boschloo test differ, two-sided'),
    0, sum((bo.pf2 <= 0.114) != bo_cells(bo.b2)), 0, "tsmethod = 'central'")

# Mehrotra, Chan and Berger (2003), Section 3 --------------------------------------------
# Two-sided tests of equality at the level 0.05, with gamma = 0.001 in the Berger-Boos
# versions (marked *). F is the Fisher test, B the Boschloo test and ZP the Z-pooled exact
# unconditional test. D and ZU are the exact unconditional tests ordered by the difference
# in proportions and by the Z statistic with the unpooled standard error. They are not
# tests of the package, since the article does not recommend them, and their p-values are
# computed here with the internal function unconditional_pvalue(). The asymptotic tests
# are the chi-squared test and the Wald test, which is the Blackwelder test with no
# margin. The p-values of the unconditional tests are refined between the grid points
# (ref.pvalue = TRUE). The null hypothesis is rejected when the p-value is below the
# level, as in BinaryRR(); the article rejects when it is at most the level, and an INFO
# item counts the p-values within 1e-6 of 0.05, for which the two rules could differ.
# The text gives the two-sided Fisher p-value as formula (2), the blaker convention. In
# the example of Section 3.1 blaker and minlike give the same p-values, and central does
# not reproduce them. The two conventions also coincide when the groups are of equal
# size, and wherever they give different rounded values in the unbalanced configurations
# of the tables, the published values of F, B and B* are those of minlike. The tables are
# therefore compared with minlike, and the note of an item gives the value of blaker when
# its rounded value differs. The percentages of Tables 1 and 3 agree with
# rounding to two decimals followed by rounding to one decimal (rule 'round2'), so that,
# for example, 2.847 is given as 2.9 and 4.945 as 5.0. An INFO item counts the items that
# ordinary rounding would reproduce and the rule 'round2' does not.
mcb.src <- 'Mehrotra, Chan and Berger (2003)'
mcb.gam <- 0.001
mcb.cols <- c('F', 'D', 'D*', 'B', 'B*', 'ZP', 'ZP*', 'ZU', 'ZU*', 'ZP asymptotic',
              'ZU asymptotic')
mcb_pv <- function(N1, N2, Test, ts, gam = 0) {
  bbssr:::get_pvalue(N1, N2, Test, 'two.sided', ts, 100L, gam, TRUE, 0)
}
# p-values of F, B and B* under the two-sided convention ts
mcb_fb <- function(N1, N2, ts) {
  list(F = mcb_pv(N1, N2, 'Fisher', ts), B = mcb_pv(N1, N2, 'Boschloo', ts),
       'B*' = mcb_pv(N1, N2, 'Boschloo', ts, mcb.gam))
}
# p-values of the eleven tests
mcb_pvalues <- function(N1, N2, ts) {
  uc <- function(stat, gam) {
    bbssr:::unconditional_pvalue(stat, N1, N2, 100L, gam, decreasing = TRUE,
                                 ref.pvalue = TRUE)
  }
  d <- abs(outer((0:N1) / N1, (0:N2) / N2, '-'))
  zu <- abs(bbssr:::zstat_margin(N1, N2, 0, 'unpooled'))
  # The outcomes with one proportion 0 and the other 1 have an infinite statistic. They
  # are the most extreme outcomes, and a finite value keeps them apart in the tie groups
  zu[is.infinite(zu)] <- max(zu[is.finite(zu)]) + 1
  fb <- mcb_fb(N1, N2, ts)
  out <- list(fb$F, uc(d, 0), uc(d, mcb.gam), fb$B, fb$`B*`,
              mcb_pv(N1, N2, 'Z-pool', ts), mcb_pv(N1, N2, 'Z-pool', ts, mcb.gam),
              uc(zu, 0), uc(zu, mcb.gam), mcb_pv(N1, N2, 'Chisq', ts),
              mcb_pv(N1, N2, 'Blackwelder', ts))
  names(out) <- mcb.cols
  out
}
# Rejection rule of BinaryRR(): p %<<% alpha of fpCompare
mcb_reject <- function(p) (0.05 - p) > sqrt(.Machine$double.eps)
# Rejection probability of a region at the response probabilities th1 and th2
mcb_rates <- function(rr, N1, N2, th1, th2 = th1) {
  b1 <- outer(0:N1, th1, function(x, t) dbinom(x, N1, t))
  b2 <- outer(0:N2, th2, function(x, t) dbinom(x, N2, t))
  colSums(b1 * (rr %*% b2))
}
# Largest rejection probability under the null hypothesis, over a grid of 4001 points with
# the three largest local maxima refined by optimize()
mcb_size <- function(rr, N1, N2) {
  g <- seq(0, 1, length.out = 4001)
  v <- mcb_rates(rr, N1, N2, g)
  peak <- which(diff(sign(diff(c(-Inf, v, -Inf)))) < 0)
  top <- head(peak[order(v[peak], decreasing = TRUE)], 3)
  ref <- vapply(top, function(k) {
    optimize(function(t) mcb_rates(rr, N1, N2, t),
             c(g[max(1, k - 1)], g[min(4001, k + 1)]), maximum = TRUE,
             tol = 1e-10)$objective
  }, numeric(1))
  max(v, ref)
}
# Rounding to two decimals followed by rounding to one decimal
mcb_round2 <- function(x) ((floor(x * 100 + 0.5 + 1e-9) + 5) %/% 10) / 10
mcb_round1 <- function(x) floor(x * 10 + 0.5 + 1e-9) / 10
# Published values of Tables 1 and 3. Table 1 lists, for each configuration (N1, N2), the
# rows theta = 0.02, 0.10, 0.25, 0.50 and the size, each with the eleven columns of
# mcb.cols. Table 3 lists theta1, theta2 and the powers of the nine exact tests
mcb.t1 <- list(
  '10, 10' = c(
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1,
    0.1, 0.1, 0.1, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 5.0,
    1.0, 1.8, 1.8, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 9.4,
    1.3, 4.1, 4.1, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 8.8,
    1.3, 4.1, 4.1, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 9.5),
  '25, 25' = c(
    0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
    0.1, 0.1, 0.1, 1.9, 1.9, 3.9, 3.9, 3.9, 3.9, 4.9, 4.9,
    2.2, 1.4, 1.4, 4.1, 4.1, 4.2, 4.2, 4.2, 4.2, 5.4, 6.4,
    3.3, 3.3, 3.3, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 6.5, 6.5,
    3.3, 3.3, 3.3, 4.6, 4.6, 4.6, 4.6, 4.6, 4.6, 6.5, 6.5),
  '50, 50' = c(
    0.0, 0.0, 0.0, 0.2, 0.2, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3,
    1.8, 0.1, 0.3, 3.2, 3.6, 3.8, 4.8, 3.8, 4.8, 5.1, 5.9,
    3.0, 1.5, 1.6, 4.1, 4.1, 4.5, 4.7, 4.5, 4.7, 5.1, 5.7,
    3.5, 3.5, 3.5, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 5.7, 5.7,
    3.5, 3.5, 3.5, 4.9, 4.9, 4.9, 4.9, 4.9, 4.9, 5.7, 6.1),
  '150, 150' = c(
    1.2, 0.0, 0.0, 2.3, 2.9, 4.6, 4.6, 4.6, 4.6, 4.6, 4.6,
    3.3, 0.1, 1.0, 3.9, 4.5, 4.8, 4.8, 4.8, 4.8, 5.0, 5.2,
    3.7, 2.0, 2.7, 4.6, 4.8, 4.8, 4.8, 4.8, 4.8, 5.0, 5.2,
    4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 5.7, 5.7,
    4.3, 4.3, 4.3, 5.0, 4.9, 5.0, 4.9, 5.0, 4.9, 5.7, 5.7),
  '16, 4' = c(
    0.2, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 5.7, 0.2,
    1.2, 0.3, 0.3, 2.9, 2.9, 2.9, 2.9, 0.0, 0.0, 8.3, 5.8,
    1.5, 1.3, 1.3, 3.7, 3.7, 3.7, 3.7, 0.4, 0.4, 4.3, 22.4,
    1.4, 3.0, 3.0, 3.4, 3.4, 3.4, 3.4, 2.8, 2.8, 5.6, 14.3,
    1.5, 3.0, 3.0, 3.9, 3.9, 3.9, 3.9, 2.8, 2.8, 9.0, 22.8),
  '40, 10' = c(
    0.8, 0.0, 0.0, 0.8, 0.8, 1.3, 1.3, 0.0, 0.0, 8.8, 0.7,
    2.4, 0.2, 0.4, 2.6, 2.6, 3.9, 3.9, 0.6, 0.6, 4.5, 20.8,
    3.4, 1.6, 1.6, 4.3, 4.3, 4.1, 4.1, 4.1, 4.1, 4.7, 9.5,
    2.9, 3.9, 3.9, 3.5, 3.5, 4.1, 4.6, 1.5, 1.5, 5.5, 8.9,
    3.4, 3.9, 3.9, 4.8, 4.8, 4.3, 4.7, 4.5, 4.5, 9.1, 21.2),
  '80, 20' = c(
    1.5, 0.0, 0.0, 1.6, 1.6, 3.3, 3.3, 0.0, 0.0, 8.7, 5.2,
    2.7, 3.4, 0.2, 3.7, 3.7, 3.2, 3.2, 3.4, 3.4, 3.7, 13.0,
    3.3, 1.7, 2.1, 4.6, 4.3, 4.6, 4.8, 1.2, 2.0, 5.1, 6.9,
    4.3, 4.0, 4.0, 4.3, 4.3, 4.1, 4.6, 0.6, 4.1, 5.1, 6.6,
    4.3, 4.0, 4.0, 5.0, 4.7, 4.6, 4.9, 3.9, 4.1, 9.2, 21.3),
  '240, 60' = c(
    1.9, 0.0, 0.2, 2.4, 3.2, 4.0, 4.0, 0.7, 0.7, 4.3, 21.3,
    3.5, 0.1, 1.2, 3.9, 4.9, 4.4, 4.7, 1.2, 2.5, 4.7, 7.0,
    4.1, 2.1, 2.9, 4.6, 4.7, 4.4, 4.8, 0.5, 4.6, 4.9, 5.7,
    4.3, 4.6, 4.6, 4.3, 4.3, 4.6, 4.6, 0.4, 4.9, 5.3, 5.5,
    4.3, 4.6, 4.6, 4.9, 4.9, 4.6, 4.9, 4.3, 4.9, 9.2, 21.4)
)
mcb.t3 <- list(
  '10, 10' = rbind(
    c(0.02, 0.54, 62.8, 66.9, 66.9, 80.8, 80.8, 80.8, 80.8, 80.8, 80.8),
    c(0.10, 0.68, 62.9, 77.7, 77.7, 79.4, 79.4, 79.4, 79.4, 79.4, 79.4),
    c(0.25, 0.84, 61.9, 78.9, 78.9, 79.2, 79.2, 79.2, 79.2, 79.2, 79.2),
    c(0.50, 0.99, 57.9, 59.9, 59.9, 78.4, 78.4, 78.4, 78.4, 78.4, 78.4)),
  '25, 25' = rbind(
    c(0.02, 0.29, 68.1, 36.8, 36.8, 76.4, 76.4, 80.4, 80.4, 80.4, 80.4),
    c(0.10, 0.45, 72.3, 66.9, 66.9, 80.9, 80.9, 80.9, 80.9, 80.9, 80.9),
    c(0.25, 0.65, 78.4, 78.4, 78.4, 80.4, 80.4, 80.4, 80.4, 80.4, 80.4),
    c(0.50, 0.86, 70.8, 69.3, 69.3, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7)),
  '50, 50' = rbind(
    c(0.02, 0.18, 68.6, 19.1, 39.6, 77.8, 78.6, 79.5, 82.2, 79.5, 82.2),
    c(0.10, 0.33, 75.7, 60.0, 65.3, 80.1, 80.3, 80.7, 82.1, 80.7, 82.1),
    c(0.25, 0.52, 74.5, 74.1, 74.1, 80.1, 80.1, 80.1, 80.1, 80.1, 80.1),
    c(0.50, 0.77, 75.3, 74.3, 74.3, 80.7, 80.7, 80.8, 80.8, 80.8, 80.8)),
  '150, 150' = rbind(
    c(0.02, 0.09, 70.4, 4.0, 43.1, 75.0, 77.7, 79.4, 79.4, 79.4, 79.4),
    c(0.10, 0.22, 77.6, 53.1, 69.0, 80.6, 81.5, 81.5, 81.5, 81.5, 81.5),
    c(0.25, 0.40, 76.1, 73.4, 73.8, 78.9, 78.5, 78.9, 78.9, 78.9, 78.9),
    c(0.50, 0.66, 78.0, 78.0, 78.0, 80.0, 79.7, 80.0, 79.7, 80.0, 79.7)),
  '16, 4' = rbind(
    c(0.02, 0.54, 64.2, 37.4, 37.4, 73.0, 73.0, 73.0, 73.0, 8.5, 8.5),
    c(0.10, 0.68, 58.3, 53.1, 53.1, 73.5, 73.5, 73.5, 73.5, 21.4, 21.4),
    c(0.25, 0.84, 47.9, 53.3, 53.3, 61.9, 61.9, 61.9, 61.9, 45.8, 45.8),
    c(0.50, 0.99, 10.1, 21.8, 21.8, 21.9, 21.9, 21.9, 21.9, 21.8, 21.8)),
  '4, 16' = rbind(
    c(0.02, 0.54, 16.3, 31.0, 31.0, 31.2, 31.2, 31.2, 31.2, 31.0, 31.0),
    c(0.10, 0.68, 41.0, 52.9, 52.9, 56.6, 56.6, 56.6, 56.6, 50.8, 50.8),
    c(0.25, 0.84, 53.7, 53.2, 53.2, 68.4, 68.4, 68.4, 68.4, 31.4, 31.4),
    c(0.50, 0.99, 63.2, 31.2, 31.2, 68.3, 68.3, 68.3, 68.3, 6.3, 6.3)),
  '40, 10' = rbind(
    c(0.02, 0.29, 68.7, 28.8, 31.5, 68.9, 68.9, 77.6, 77.6, 4.0, 4.0),
    c(0.10, 0.45, 67.3, 46.6, 50.0, 71.3, 71.3, 71.8, 71.8, 16.6, 16.6),
    c(0.25, 0.65, 60.0, 59.8, 59.8, 68.7, 68.7, 65.7, 66.9, 29.1, 29.1),
    c(0.50, 0.86, 50.2, 51.7, 51.7, 54.0, 54.0, 56.1, 56.4, 42.9, 42.9)),
  '10, 40' = rbind(
    c(0.02, 0.29, 30.3, 12.9, 12.9, 30.5, 30.5, 30.5, 30.5, 70.4, 70.4),
    c(0.10, 0.45, 51.2, 48.6, 48.6, 56.1, 56.1, 56.8, 56.8, 47.1, 47.1),
    c(0.25, 0.65, 55.6, 60.9, 60.9, 60.9, 60.9, 61.9, 65.2, 36.8, 36.8),
    c(0.50, 0.86, 64.1, 49.5, 50.1, 68.7, 68.7, 68.5, 68.5, 18.5, 18.5)),
  '80, 20' = rbind(
    c(0.02, 0.18, 64.6, 12.9, 38.5, 70.6, 70.6, 76.2, 76.2, 2.4, 2.4),
    c(0.10, 0.33, 65.0, 39.9, 51.1, 68.1, 68.1, 68.4, 68.4, 10.6, 10.6),
    c(0.25, 0.52, 57.4, 54.6, 55.7, 64.8, 63.3, 63.9, 63.9, 18.4, 41.8),
    c(0.50, 0.77, 59.5, 56.2, 56.2, 59.7, 59.7, 57.7, 59.9, 32.1, 60.7)),
  '20, 80' = rbind(
    c(0.02, 0.18, 32.4, 2.9, 5.0, 40.6, 39.8, 40.7, 40.7, 62.7, 62.7),
    c(0.10, 0.33, 49.4, 39.5, 42.1, 56.5, 55.7, 53.2, 56.9, 42.6, 57.1),
    c(0.25, 0.52, 58.6, 56.0, 56.0, 58.6, 58.6, 56.9, 59.0, 30.1, 59.0),
    c(0.50, 0.77, 59.3, 54.5, 56.3, 65.9, 64.8, 65.7, 65.7, 18.2, 37.8)),
  '240, 60' = rbind(
    c(0.02, 0.09, 64.1, 3.4, 38.2, 67.3, 68.5, 68.9, 69.2, 2.8, 2.8),
    c(0.10, 0.22, 62.7, 33.2, 53.4, 65.1, 67.4, 67.8, 67.8, 13.4, 41.4),
    c(0.25, 0.40, 59.8, 53.4, 56.2, 61.1, 61.2, 61.1, 62.3, 18.1, 54.5),
    c(0.50, 0.66, 59.2, 59.6, 59.6, 59.2, 59.2, 59.7, 60.0, 25.3, 61.5)),
  '60, 240' = rbind(
    c(0.02, 0.09, 41.1, 0.1, 8.5, 41.1, 48.5, 43.2, 44.6, 45.0, 46.8),
    c(0.10, 0.22, 54.1, 31.5, 43.5, 57.5, 57.8, 55.5, 58.3, 37.0, 66.2),
    c(0.25, 0.40, 56.0, 54.5, 54.6, 58.6, 58.6, 57.5, 58.5, 27.1, 61.8),
    c(0.50, 0.66, 58.7, 59.0, 59.0, 62.7, 62.0, 61.3, 62.2, 21.0, 59.0))
)
# Section 3.1: 8 of 148 against 1 of 132 responders. Two-sided p-values and the exact
# 99.9 per cent Clopper-Pearson interval for the common response probability
mcb.ex <- mcb_pvalues(148, 132, 'blaker')
mcb.lab <- 'MCB2003 Section 3.1: 8 / 148 against 1 / 132'
add(mcb.src, paste0(mcb.lab, ': 99.9% confidence interval, ', c('lower', 'upper')),
    c(0.0080, 0.0826), unname(bbssr:::cp_bounds(280, mcb.gam)[10, ]), 4)
# The p-values of F, B and B* under minlike, and of F under central, for the notes
mcb.exm <- vapply(mcb_fb(148, 132, 'minlike'), function(p) p[9, 2], numeric(1))
mcb.exc <- mcb_pv(148, 132, 'Fisher', 'central')[9, 2]
mcb.exnote <- c(sprintf("tsmethod = 'blaker'; minlike gives %.6f and central %.6f",
                        mcb.exm[['F']], mcb.exc), '', '',
                sprintf("tsmethod = 'blaker'; minlike gives %.6f", mcb.exm[c('B', 'B*')]),
                rep('', 4))
add(mcb.src, paste0(mcb.lab, ': two-sided p-value, ', mcb.cols[1:9]),
    c(0.0388, 0.4386, 0.1603, 0.0347, 0.0325, 0.0291, 0.0282, 0.0229, 0.0215),
    vapply(mcb.ex[1:9], function(p) p[9, 2], numeric(1)), 4, mcb.exnote)
# Tables 1 (type I error rates and sizes) and 3 (powers), in per cent
mcb.rows1 <- c(sprintf('theta = %.2f', c(0.02, 0.10, 0.25, 0.50)), 'size')
mcb.fbi <- match(c('F', 'B', 'B*'), mcb.cols)
mcb.near <- 0
mcb.tab <- list()
for (key in union(names(mcb.t1), names(mcb.t3))) {
  n <- as.integer(strsplit(key, ', ')[[1]])
  N1 <- n[1]
  N2 <- n[2]
  P <- mcb_pvalues(N1, N2, 'minlike')
  mcb.near <- mcb.near + sum(vapply(P, function(p) {
    as.numeric(sum(abs(p - 0.05) < 1e-6))
  }, numeric(1)))
  RR <- lapply(P, mcb_reject)
  # The blaker convention differs from minlike only when the groups are of unequal size
  RRb <- if (N1 != N2) lapply(mcb_fb(N1, N2, 'blaker'), mcb_reject) else RR[mcb.fbi]
  note_fb <- function(rec, recb) {
    nt <- matrix('', nrow(rec), ncol(rec))
    nt[, mcb.fbi] <- ifelse(mcb_round2(recb) != mcb_round2(rec[, mcb.fbi, drop = FALSE]),
                            sprintf('blaker (formula (2)) gives %.3f', recb), '')
    nt
  }
  if (!is.null(mcb.t1[[key]])) {
    th <- c(0.02, 0.10, 0.25, 0.50)
    rate <- function(rr) c(mcb_rates(rr, N1, N2, th), mcb_size(rr, N1, N2)) * 100
    rec <- vapply(RR, rate, numeric(5))
    recb <- vapply(RRb, rate, numeric(5))
    mcb.tab[[length(mcb.tab) + 1L]] <- data.frame(
      item = c(outer(mcb.rows1, mcb.cols, function(r, k) {
        sprintf('MCB2003 Table 1: (%d, %d), %s: %s', N1, N2, r, k)
      })),
      published = c(matrix(mcb.t1[[key]], nrow = 5, byrow = TRUE)),
      recomputed = c(rec), note = c(note_fb(rec, recb)), stringsAsFactors = FALSE)
  }
  if (!is.null(mcb.t3[[key]])) {
    t3 <- mcb.t3[[key]]
    pw <- function(rr) {
      vapply(seq_len(nrow(t3)), function(i) {
        mcb_rates(rr, N1, N2, t3[i, 1], t3[i, 2])
      }, numeric(1)) * 100
    }
    rec <- matrix(vapply(RR[1:9], pw, numeric(nrow(t3))), nrow = nrow(t3))
    recb <- matrix(vapply(RRb, pw, numeric(nrow(t3))), nrow = nrow(t3))
    rows3 <- sprintf('(%.2f, %.2f)', t3[, 1], t3[, 2])
    mcb.tab[[length(mcb.tab) + 1L]] <- data.frame(
      item = c(outer(rows3, mcb.cols[1:9], function(r, k) {
        sprintf('MCB2003 Table 3: (%d, %d), %s: %s', N1, N2, r, k)
      })),
      published = c(t3[, 3:11]), recomputed = c(rec), note = c(note_fb(rec, recb)),
      stringsAsFactors = FALSE)
  }
}
mcb.tab <- do.call(rbind, mcb.tab)
mcb.miss <- abs(mcb_round2(mcb.tab$recomputed) - mcb.tab$published) > 1e-9
mcb.tab$note[mcb.miss] <- trimws(paste(
  'Not reproduced, cause not identified. The independent implementation in',
  'tools/reference/check_mehrotra_2003.py gives the same value.', mcb.tab$note[mcb.miss]))
add(mcb.src, mcb.tab$item, mcb.tab$published, mcb.tab$recomputed, 1, mcb.tab$note,
    rule = 'round2')
add(mcb.src, 'MCB2003 Tables 1 and 3: p-values within 1e-6 of the level 0.05', NA,
    mcb.near, NA, paste('INFO: count over the p-values used for the tables, with minlike',
                        'for F, B and B*'))
add(mcb.src, paste('MCB2003 Tables 1 and 3: items that ordinary rounding reproduces and',
                   'the rule round2 does not'), NA,
    sum(abs(mcb_round1(mcb.tab$recomputed) - mcb.tab$published) < 1e-9 & mcb.miss), NA,
    'INFO: check of the rounding rule')

# Berger and Boos (1994), Example 2 ------------------------------------------------------
# 14 of 47 against 48 of 283 responders. The p-value of the two-sided Z-pooled test, whose
# ordering is that of the Pearson chi-squared statistic, maximized over [0, 1] and over
# the .999 confidence interval, to which gamma = .001 is added. The statistic and the
# interval agree with the article after rounding. The maximum over [0, 1] (.061), the
# maximum over the interval (.036) and p_.001 (.037) are compared after truncation to
# three decimals: rounding gives .037 and .038 for the last two, and the article also
# writes the location 0.0039 of the maximum over [0, 1] as p(.003). The tail probability
# is symmetric about 1/2, so the maximum is attained at that point and at its mirror
# image; the location in [0, 1/2] is reported as INFO
bb.src <- 'Berger and Boos (1994)'
bb.lab <- 'BB1994 Example 2: 14 / 47 against 48 / 283'
bb.z <- bbssr:::zstat(47, 283)
bb.trunc <- 'Truncated to three decimals, as the article does for the maximized values'
add(bb.src, paste0(bb.lab, ': chi-squared statistic'), 4.346, bb.z[15, 49] ^ 2, 3)
add(bb.src, paste0(bb.lab, ': .999 confidence interval, ', c('lower', 'upper')),
    c(0.123, 0.267), unname(bbssr:::cp_bounds(330, 0.001)[63, ]), 3)
bb.sup <- bbssr:::get_pvalue(47, 283, 'Z-pool', 'two.sided', 'minlike', 100L, 0, TRUE,
                             0)[15, 49]
bb.p <- bbssr:::get_pvalue(47, 283, 'Z-pool', 'two.sided', 'minlike', 100L, 0.001, TRUE,
                           0)[15, 49]
add(bb.src, paste0(bb.lab, c(': p-value maximized over [0, 1]',
                             ': maximum over the confidence interval',
                             ': p-value p_.001')),
    c(0.061, 0.036, 0.037), c(bb.sup, bb.p - 0.001, bb.p), 3, bb.trunc,
    rule = 'truncate')
bb.mask <- abs(bb.z) >= abs(bb.z[15, 49]) * (1 - 1e-10)
bb.tail <- function(t) mcb_rates(bb.mask, 47, 283, t)
bb.th <- seq(0, 0.5, by = 1e-4)
bb.k <- which.max(bb.tail(bb.th))
bb.loc <- optimize(bb.tail, c(bb.th[max(1, bb.k - 1)], bb.th[min(5001, bb.k + 1)]),
                   maximum = TRUE, tol = 1e-10)$maximum
add(bb.src, paste0(bb.lab, ': location of the maximum in [0, 1/2]'), 0.003, bb.loc, NA,
    paste('INFO: grid of step 1e-4 refined by optimize(); the maximum is attained again',
          'at the mirror image about 1/2, and the article gives p(.003) = .061'))

# Fay and Hunsberger (2021), Section 8 and Table 1 ---------------------------------------
# 8 of 14 against 1 of 7 responders. Two-sided Fisher p-values under the three
# conventions, and the ordering function of Blaker's test, T_B(x, 1), which is the blaker
# p-value of each outcome with 9 responders in total
fh.src <- 'Fay and Hunsberger (2021)'
fh.pub <- c(blaker = 0.087, minlike = 0.159, central = 0.157)
fh.p <- lapply(names(fh.pub), function(ts) {
  bbssr:::get_pvalue(14, 7, 'Fisher', 'two.sided', ts, 100L, 0, FALSE, 0)
})
add(fh.src, paste0('FH2021 Section 8: 8 / 14 against 1 / 7: two-sided p-value, ',
                   names(fh.pub)), fh.pub, vapply(fh.p, function(p) p[9, 2], numeric(1)),
    3)
add(fh.src, sprintf('FH2021 Table 1: T_B(x, 1) at x2 = %d', 0:7),
    c(0.007, 0.087, 0.642, 1.000, 0.397, 0.159, 0.016, 0.000),
    fh.p[[1]][cbind(9 - 0:7 + 1, 0:7 + 1)], 3)

# Verdicts and summary -------------------------------------------------------------------
tab <- do.call(rbind, cmp)
info <- is.na(tab$digits) | is.na(tab$published)
d <- ifelse(info, 0, tab$digits)
# The rule 'round2' rounds to one more decimal first and then to the digits of the
# publication
shown <- ifelse(tab$rule == 'truncate', floor(tab$recomputed * 10^d + 1e-9) / 10^d,
                ifelse(tab$rule == 'round2',
                       ((floor(tab$recomputed * 10^(d + 1) + 0.5 + 1e-9) + 5) %/% 10) /
                         10^d,
                       round(tab$recomputed, d)))
tab$verdict <- ifelse(info, 'INFO',
                      ifelse(abs(shown - tab$published) < 1e-9, 'PASS', 'FAIL'))
tab$tolerance <- NA_real_
for (e in explained) {
  hit <- tab$verdict == 'FAIL' & tab$item %in% e$items &
    abs(tab$recomputed - tab$published) <= e$tol + 1e-12
  tab$verdict[hit] <- 'EXPLAINED'
  tab$tolerance[hit] <- e$tol
  tab$note[hit] <- e$reason
}
write.csv(tab, file.path(out.dir, 'published-comparison.csv'), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out.dir, 'session.txt'))

counts <- table(factor(tab$verdict, levels = c('PASS', 'EXPLAINED', 'FAIL', 'INFO')))
fmt <- function(x) ifelse(is.na(x), '', formatC(x, digits = 6, format = 'g'))
md <- c(
  '# Reproduction of published results with bbssr',
  '',
  sprintf('Generated by inst/reproduce/reproduce-published.R on %s (bbssr %s).',
          format(Sys.time(), '%Y-%m-%d %H:%M'), as.character(utils::packageVersion('bbssr'))),
  sprintf('Run time: %.1f minutes.', as.numeric(difftime(Sys.time(), t.start,
                                                         units = 'mins'))),
  '',
  sprintf('PASS: %d, EXPLAINED: %d, FAIL: %d, INFO: %d', counts[['PASS']],
          counts[['EXPLAINED']], counts[['FAIL']], counts[['INFO']]),
  '',
  'A value passes when the recomputed value, rounded to the digits of the publication',
  '(truncated where the publication truncates, and rounded first to one more digit where',
  'the publication does so), equals the published value. A value that',
  'does not pass is EXPLAINED when a documented reason covers it and its difference from',
  'the published value is within the tolerance stated with the reason; otherwise it is a',
  'FAIL.',
  '',
  '| Verdict | Item | Published | Recomputed | Tolerance | Note |',
  '|---|---|---|---|---|---|',
  sprintf('| %s | %s | %s | %s | %s | %s |', tab$verdict, tab$item, fmt(tab$published),
          fmt(tab$recomputed), fmt(tab$tolerance), tab$note)
)
writeLines(md, file.path(out.dir, 'summary.md'))
print(counts)
# Item names must be unique, and every explained item must have a comparison. This is
# checked after the output is written, so that a mistake does not discard the results
if (anyDuplicated(tab$item)) {
  dup <- unique(tab$item[duplicated(tab$item)])
  stop('duplicated items: ', paste(dup, collapse = '; '))
}
covered <- unlist(lapply(explained, function(e) e$items))
if (!all(covered %in% tab$item)) {
  stop('explained items without a comparison: ',
       paste(setdiff(covered, tab$item), collapse = '; '))
}
