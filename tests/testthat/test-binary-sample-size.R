test_that("BinarySampleSize returns a bbssr_samplesize data frame", {
  res <- BinarySampleSize(p1 = 0.5, p2 = 0.2, r = 1, alpha = 0.025,
                          tar.power = 0.8, Test = 'Chisq')
  expect_s3_class(res, 'bbssr_samplesize')
  expect_equal(nrow(res), 1L)
  expect_named(res, c('p1', 'p2', 'r', 'alpha', 'tar.power', 'Test', 'alternative',
                      'Power', 'N1', 'N2', 'N'))
  expect_type(res$N1, 'integer')
  expect_type(res$N2, 'integer')
  expect_equal(res$N, res$N1 + res$N2)
})

test_that("the returned sample size is the smallest one attaining the target power", {
  for (tst in c('Chisq', 'Fisher')) {
    res <- BinarySampleSize(0.6, 0.25, 1, 0.025, 0.8, tst)
    at <- BinaryPower(0.6, 0.25, res$N1, res$N2, 0.025, tst)$Power
    below <- BinaryPower(0.6, 0.25, ceiling(res$N2 - 1), res$N2 - 1, 0.025, tst)$Power
    expect_gte(at, 0.8)
    expect_lt(below, 0.8)
  }
})

test_that("the allocation ratio is respected", {
  for (r in c(0.5, 1, 2)) {
    res <- BinarySampleSize(0.6, 0.3, r, 0.025, 0.8, 'Chisq')
    expect_equal(res$N1, as.integer(ceiling(r * res$N2)))
  }
})

test_that("a larger target power requires a larger sample size", {
  n <- vapply(c(0.7, 0.8, 0.9), function(tp) {
    BinarySampleSize(0.5, 0.25, 1, 0.025, tp, 'Chisq')$N
  }, numeric(1))
  expect_true(all(diff(n) > 0))
})

test_that("a smaller treatment effect requires a larger sample size", {
  n <- vapply(c(0.7, 0.6, 0.5), function(p1) {
    BinarySampleSize(p1, 0.3, 1, 0.025, 0.8, 'Chisq')$N
  }, numeric(1))
  expect_true(all(diff(n) > 0))
})

test_that("the two-sided test requires at least as many patients as the one-sided test", {
  one <- BinarySampleSize(0.6, 0.3, 1, 0.025, 0.8, 'Fisher')$N
  two <- BinarySampleSize(0.6, 0.3, 1, 0.025, 0.8, 'Fisher',
                          alternative = 'two.sided')$N
  expect_gte(two, one)
})

test_that("the Boschloo test needs no more patients than the Fisher test", {
  fisher <- BinarySampleSize(0.6, 0.2, 1, 0.025, 0.8, 'Fisher')$N
  boschloo <- BinarySampleSize(0.6, 0.2, 1, 0.025, 0.8, 'Boschloo')$N
  expect_lte(boschloo, fisher)
})

test_that("BinarySampleSize validates its arguments", {
  expect_error(BinarySampleSize(0.4, 0.4, 1, 0.025, 0.8, 'Chisq'), 'must differ')
  expect_error(BinarySampleSize(0.5, 0.2, 0, 0.025, 0.8, 'Chisq'), 'positive')
  expect_error(BinarySampleSize(0.5, 0.2, 1, 0.025, 1, 'Chisq'), 'tar.power')
  expect_error(BinarySampleSize(c(0.5, 0.6), 0.2, 1, 0.025, 0.8, 'Chisq'), 'single value')
})

test_that("the normal approximation gives the published sample sizes", {
  # Kieser (2020), Example 21.1: 2 x 157 patients
  k <- BinarySampleSize(0.42, 0.27, 1, 0.05, 0.8, 'Chisq', alternative = 'two.sided',
                        method = 'standard')
  expect_equal(c(k$N1, k$N2), c(157L, 157L))
  # Friede and Kieser (2004), Table I. The larger group receives the larger probability
  fk <- function(p, r, rounding) {
    BinarySampleSize(p + 0.2, p, r, 0.05, 0.8, 'Chisq', alternative = 'two.sided',
                     method = 'null.variance', rounding = rounding)$N
  }
  expect_equal(vapply(c(0.05, 0.2, 0.4), fk, numeric(1), r = 1, rounding = 'group'),
               c(102, 166, 198))
  expect_equal(vapply(c(0.05, 0.2, 0.4, 0.6, 0.75), fk, numeric(1), r = 3,
                      rounding = 'friede-kieser'),
               c(168, 239, 260, 198, 95))
})

test_that("the lower-tail alternative mirrors the upper-tail one", {
  for (tst in c('Chisq', 'Boschloo')) {
    up <- BinarySampleSize(0.6, 0.3, 1, 0.025, 0.8, tst)
    down <- BinarySampleSize(0.3, 0.6, 1, 0.025, 0.8, tst, alternative = 'less')
    expect_equal(c(down$N1, down$N2), c(up$N1, up$N2), info = tst)
    expect_equal(down$Power, up$Power, tolerance = 1e-12, info = tst)
  }
})

test_that("the direction of the effect must agree with the alternative", {
  expect_error(BinarySampleSize(0.2, 0.5, 1, 0.025, 0.8, 'Chisq'), "'greater'")
  expect_error(BinarySampleSize(0.5, 0.2, 1, 0.025, 0.8, 'Chisq', alternative = 'less'),
               "'less'")
  expect_error(BinarySampleSize(0.5, 0.2, 1, 0.025, 0.8, 'Chisq', rounding = 'total'),
               'rounding')
})

test_that("BinarySampleSize passes ref.pvalue to every p-value computation", {
  seen <- ref_pvalue_calls(BinarySampleSize(0.6, 0.25, 1, 0.025, 0.8, 'Z-pool',
                                            ref.pvalue = TRUE))
  expect_gt(nrow(seen), 1)
  expect_true(all(seen$ref.pvalue))
  ss <- BinarySampleSize(0.6, 0.25, 1, 0.025, 0.8, 'Z-pool', ref.pvalue = TRUE)
  expect_true(attr(ss, 'ref.pvalue'))
  seen <- ref_pvalue_calls(BinarySampleSize(0.6, 0.25, 1, 0.025, 0.8, 'Z-pool'))
  expect_gt(nrow(seen), 0)
  expect_false(any(seen$ref.pvalue))
})

test_that("the sample sizes of Farrington and Manning (1990, Table I) are reproduced", {
  # One-sided level 0.05, power 0.9 and each group rounded to the nearest whole number.
  # The article writes theta = N2 / N1, so r = 1 / theta, and tests s <= s0, so the
  # margin is -s0
  tab <- data.frame(p1 = c(0.1, 0.2, 0.5, 0.1, 0.25), p2 = c(0.1, 0.1, 0.1, 0.05, 0.05),
                    s0 = c(-0.2, -0.1, 0.2, -0.05, 0.1),
                    theta = c(2 / 3, 1, 3 / 2, 2 / 3, 1),
                    N1 = c(63, 57, 67, 168, 197), N2 = c(42, 57, 101, 112, 197))
  for (i in seq_len(nrow(tab))) {
    ss <- BinarySampleSize(tab$p1[i], tab$p2[i], 1 / tab$theta[i], 0.05, 0.9,
                           'Farrington-Manning', method = 'standard',
                           rounding = 'nearest', margin = -tab$s0[i])
    expect_equal(c(ss$N1, ss$N2), c(tab$N1[i], tab$N2[i]), info = i)
  }
})

test_that("the fixed sample sizes of Friede et al. (2007, Table 2) are reproduced", {
  # Farrington-Manning test, one-sided level 0.025, power 0.8, margin 0.1, p1 = p2 and the
  # two groups rounded up separately. The article writes q = n2 / n1, so r = 1 / q
  pa <- c(0.3, 0.4, 0.5, 0.6, 0.7)
  q <- c(1 / 3, 1 / 2, 1)
  want <- rbind(c(928, 770, 658), c(1023, 857, 750), c(1035, 876, 780), c(966, 825, 750),
                c(815, 707, 658))
  for (i in seq_along(pa)) {
    for (j in seq_along(q)) {
      N <- BinarySampleSize(pa[i], pa[i], 1 / q[j], 0.025, 0.8, 'Farrington-Manning',
                            method = 'standard', rounding = 'friede-kieser',
                            margin = 0.1)$N
      expect_equal(N, want[i, j], info = sprintf('pa = %g, q = %g', pa[i], q[j]))
    }
  }
})

test_that("the sample sizes of Blackwelder (1982, Table 3) are reproduced", {
  # One-sided level 0.05, power 0.9, equal groups and the size of each group rounded up.
  # Group 1 is the experimental therapy. The entry 0.6, 0.55, 0.1 is left out: the article
  # gives 3342 from the rounded quantiles 1.645 and 1.282, and the exact quantiles give
  # 3340
  ni <- data.frame(ps = c(0.9, 0.9, 0.6, 0.6, 0.4, 0.4, 0.9, 0.9, 0.6),
                   pe = c(0.9, 0.9, 0.6, 0.6, 0.4, 0.4, 0.85, 0.8, 0.5),
                   delta = c(0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.2),
                   N = c(310, 78, 824, 206, 824, 206, 1492, 430, 840))
  for (i in seq_len(nrow(ni))) {
    N <- BinarySampleSize(ni$pe[i], ni$ps[i], 1, 0.05, 0.9, 'Blackwelder',
                          method = 'alternative.variance', margin = ni$delta[i])$N
    expect_equal(N, ni$N[i], info = i)
  }
  # The conventional null hypothesis, with group 1 the standard therapy
  sup <- data.frame(ps = c(0.9, 0.9, 0.6, 0.6, 0.4, 0.4),
                    pe = c(0.8, 0.7, 0.5, 0.4, 0.3, 0.2),
                    N = c(430, 130, 840, 206, 772, 172))
  for (i in seq_len(nrow(sup))) {
    N <- BinarySampleSize(sup$ps[i], sup$pe[i], 1, 0.05, 0.9, 'Blackwelder',
                          method = 'alternative.variance')$N
    expect_equal(N, sup$N[i], info = i)
  }
})

test_that("the exact search handles a non-inferiority margin", {
  # Reference values from tools/reference/reference_values.py, whose search starts from
  # the formula of Farrington and Manning and moves one patient at a time
  seen <- ref_pvalue_calls(ss <- BinarySampleSize(0.8, 0.8, 1, 0.025, 0.8,
                                                  'Farrington-Manning', margin = 0.15))
  expect_equal(ss$N2, 113)
  expect_gt(nrow(seen), 0)
  expect_true(all(seen$margin == 0.15))
  expect_equal(attr(ss, 'margin'), 0.15)
  expect_equal(BinarySampleSize(0.6, 0.65, 2, 0.025, 0.9, 'Blackwelder', margin = 0.2)$N2,
               160)
})

test_that("BinarySampleSize checks the direction of the effect against the margin", {
  expect_error(BinarySampleSize(0.6, 0.75, 1, 0.025, 0.8, 'Farrington-Manning',
                                margin = 0.1), 'exceed -margin')
  expect_error(BinarySampleSize(0.75, 0.6, 1, 0.025, 0.8, 'Farrington-Manning',
                                alternative = 'less', margin = 0.1), 'fall below margin')
  expect_error(BinarySampleSize(0.6, 0.6, 1, 0.025, 0.8, 'Chisq', margin = 0.1),
               'Blackwelder')
  expect_no_error(BinarySampleSize(0.6, 0.6, 1, 0.025, 0.8, 'Blackwelder',
                                   method = 'alternative.variance', margin = 0.1))
})
