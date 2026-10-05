# Reproduces the published numerical results of Friede and Kieser (2004, Table I and
# Section 5), Kieser (2020, Example 21.1), Farrington and Manning (1990, Tables I and II and
# the first example), Blackwelder (1982, Table 3 and the examples) and Friede, Mitchell and
# Mueller-Velten (2007, Tables 2 and 3 and Sections 5 and 6) with bbssr.
# Run from the package root after devtools::load_all(). The results are written to
# reproduce-output/: published-comparison.csv lists every published value next to the
# recomputed value with a verdict (PASS, EXPLAINED, FAIL or INFO), and summary.md is
# generated from it. The full run takes several minutes.

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
explain('Kieser 2020: fixed design, Delta = 0.30, maximum level', 1e-4, paste(
  'The range of the overall rate behind the summary statistics is not stated. Over',
  '[0.15, 0.85] the recomputed maximum is reported; outside this range the level of the',
  'fixed design rises to about 0.029 near p = 0.09.'))
explain('Kieser 2020: IPS design, Delta = 0.30, maximum level', 1e-4, paste(
  'Same as for the fixed design: the range of the overall rate is not stated, and the',
  'recomputed maximum is taken over [0.15, 0.85].'))

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
tie <- do.call(BinaryTypeIErrorBSSR, c(k.args, list(refine = FALSE)))
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
in.range <- k.grid >= 0.05 - 1e-9 & k.grid <= 0.95 + 1e-9
add('Kieser (2020)', 'Kieser 2020: fixed design, Delta = 0.15, mean level', 0.0248,
    mean(tie$TIE.TRAD[in.range]), NA,
    'INFO: range of the overall rate not stated; recomputed over [0.05, 0.95]')
add('Kieser (2020)', 'Kieser 2020: IPS design, Delta = 0.15, mean level', 0.0250,
    mean(tie$TIE.BSSR[in.range]), NA,
    'INFO: range of the overall rate not stated; recomputed over [0.05, 0.95]')
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
                                      n.interim = c(20, 20), refine = FALSE)))
write.csv(k3, file.path(out.dir, 'kieser-2020-level-delta030.csv'), row.names = FALSE)
mid <- k.grid >= 0.15 - 1e-9 & k.grid <= 0.85 + 1e-9
range.note <- 'Range of the overall rate not stated; recomputed over [0.15, 0.85]'
add('Kieser (2020)', 'Kieser 2020: fixed design, Delta = 0.30, maximum level', 0.0268,
    max(k3$TIE.TRAD[mid]), 4, range.note)
add('Kieser (2020)', 'Kieser 2020: fixed design, Delta = 0.30, minimum level', 0.0253,
    min(k3$TIE.TRAD[mid]), 4, range.note)
add('Kieser (2020)', 'Kieser 2020: IPS design, Delta = 0.30, maximum level', 0.0273,
    max(k3$TIE.BSSR[mid]), 4, range.note)
add('Kieser (2020)', 'Kieser 2020: fixed design, Delta = 0.30, maximum over [0.01, 0.99]',
    NA, max(k3$TIE.TRAD), NA, 'INFO')
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
add(fm.src, paste0(fm.lab, ': true power'), fm$pw.pub, fm$power.article, 4,
    ifelse(fm$power - fm$power.article > 5e-5,
           sprintf('Without the outcome (0, 0); bbssr, which rejects it, gives %.4f',
                   fm$power), ''),
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

# Verdicts and summary -------------------------------------------------------------------
tab <- do.call(rbind, cmp)
info <- is.na(tab$digits) | is.na(tab$published)
d <- ifelse(info, 0, tab$digits)
shown <- ifelse(tab$rule == 'truncate', floor(tab$recomputed * 10^d + 1e-9) / 10^d,
                round(tab$recomputed, d))
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
  '(truncated where the publication truncates), equals the published value. A value that',
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
