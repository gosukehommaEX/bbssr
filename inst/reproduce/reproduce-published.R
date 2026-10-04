# Reproduces the published numerical results of Friede and Kieser (2004, Table I and
# Section 5) and Kieser (2020, Example 21.1) with bbssr.
# Run from the package root after devtools::load_all(). The results are written to
# reproduce-output/: published-comparison.csv lists every published value next to the
# recomputed value with a verdict (PASS, EXPLAINED, FAIL or INFO), and summary.md is
# generated from it. The full run takes a few minutes.

out.dir <- 'reproduce-output'
dir.create(out.dir, showWarnings = FALSE)
t.start <- Sys.time()
cmp <- list()
add <- function(source, item, published, recomputed, digits, note = '') {
  cmp[[length(cmp) + 1L]] <<- data.frame(source = source, item = item,
                                         published = published, recomputed = recomputed,
                                         digits = digits, note = note,
                                         stringsAsFactors = FALSE)
}
# Discrepancies that are understood, with the reason. A row listed here is reported as
# EXPLAINED instead of FAIL; it is never used to hide an unexpected difference
explained <- list(
  'Kieser 2020: fixed design, Delta = 0.30, maximum level' = paste(
    'The range of the overall rate behind the summary statistics is not stated. Over',
    '[0.15, 0.85] the recomputed maximum is reported; outside this range the level of the',
    'fixed design rises to about 0.029 near p = 0.09.'),
  'Kieser 2020: IPS design, Delta = 0.30, maximum level' = paste(
    'Same as for the fixed design: the range of the overall rate is not stated, and the',
    'recomputed maximum is taken over [0.15, 0.85].')
)

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

# Verdicts and summary -------------------------------------------------------------------
tab <- do.call(rbind, cmp)
tab$verdict <- ifelse(
  is.na(tab$digits) | is.na(tab$published), 'INFO',
  ifelse(abs(round(tab$recomputed, tab$digits) - tab$published) < 1e-9, 'PASS',
         ifelse(tab$item %in% names(explained), 'EXPLAINED', 'FAIL'))
)
hit <- tab$verdict == 'EXPLAINED'
tab$note[hit] <- unlist(explained[tab$item[hit]])
write.csv(tab, file.path(out.dir, 'published-comparison.csv'), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out.dir, 'session.txt'))

counts <- table(factor(tab$verdict, levels = c('PASS', 'EXPLAINED', 'FAIL', 'INFO')))
fmt <- function(x) ifelse(is.na(x), '', format(signif(x, 6)))
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
  'A value passes when the recomputed value, rounded to the digits of the publication,',
  'equals the published value.',
  '',
  '| Verdict | Item | Published | Recomputed | Note |',
  '|---|---|---|---|---|',
  sprintf('| %s | %s | %s | %s | %s |', tab$verdict, tab$item, fmt(tab$published),
          fmt(tab$recomputed), tab$note)
)
writeLines(md, file.path(out.dir, 'summary.md'))
print(counts)
