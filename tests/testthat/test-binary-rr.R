all_tests <- c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')
exact_tests <- c('Fisher', 'Z-pool', 'Boschloo')

test_that("BinaryRR returns a bbssr_rr object of the right shape", {
  for (tst in all_tests) {
    RR <- BinaryRR(N1 = 4, N2 = 5, alpha = 0.05, Test = tst, n.grid = 20)
    expect_s3_class(RR, 'bbssr_rr')
    expect_true(is.matrix(RR))
    expect_type(as.vector(RR), 'logical')
    expect_equal(dim(RR), c(5L, 6L))
    expect_equal(attr(RR, 'Test'), tst)
    expect_equal(dim(attr(RR, 'p.value')), c(5L, 6L))
  }
})

test_that("BinaryRR reproduces the chi-squared rejection region analytically", {
  N1 <- 12
  N2 <- 9
  alpha <- 0.025
  RR <- BinaryRR(N1, N2, alpha, 'Chisq')
  expect_equal(as_plain(RR), zstat(N1, N2) > stats::qnorm(1 - alpha))
  RR2 <- BinaryRR(N1, N2, alpha, 'Chisq', alternative = 'two.sided')
  expect_equal(as_plain(RR2), abs(zstat(N1, N2)) > stats::qnorm(1 - alpha / 2))
})

test_that("BinaryRR rejects only in the direction of the alternative", {
  N1 <- 10
  N2 <- 10
  rd <- outer(0:N1 / N1, 0:N2 / N2, '-')
  for (tst in all_tests) {
    RR <- as_plain(BinaryRR(N1, N2, 0.05, tst, n.grid = 30))
    expect_true(all(rd[RR] > 0), info = tst)
  }
})

test_that("the one-sided rejection region is monotone in the outcome grid", {
  # More responders in group 1 or fewer in group 2 can only help rejection
  for (tst in all_tests) {
    RR <- as_plain(BinaryRR(8, 8, 0.05, tst, n.grid = 30))
    for (j in seq_len(ncol(RR))) {
      expect_true(all(diff(RR[, j]) >= 0), info = sprintf('%s, column %d', tst, j))
    }
    for (i in seq_len(nrow(RR))) {
      expect_true(all(diff(RR[i, ]) <= 0), info = sprintf('%s, row %d', tst, i))
    }
  }
})

test_that("the two-sided rejection region is symmetric when the groups are balanced", {
  N <- 8
  for (tst in all_tests) {
    for (ts in c('minlike', 'central', 'blaker')) {
      RR <- as_plain(BinaryRR(N, N, 0.05, tst, alternative = 'two.sided',
                              tsmethod = ts, n.grid = 30))
      expect_equal(RR, t(RR), info = sprintf('%s, %s', tst, ts))
    }
  }
})

test_that("the central two-sided region is the union of the two one-sided regions", {
  N1 <- 9
  N2 <- 7
  alpha <- 0.02
  for (tst in c('Fisher', 'Fisher-midP')) {
    two <- as_plain(BinaryRR(N1, N2, 2 * alpha, tst,
                             alternative = 'two.sided', tsmethod = 'central'))
    up <- as_plain(BinaryRR(N1, N2, alpha, tst))
    lo <- t(as_plain(BinaryRR(N2, N1, alpha, tst)))
    expect_equal(two, up | lo, info = tst)
  }
})

test_that("the exact tests control the type I error rate", {
  alpha <- 0.05
  for (tst in exact_tests) {
    for (alt in c('greater', 'two.sided')) {
      RR <- BinaryRR(8, 8, alpha, tst, alternative = alt, n.grid = 60)
      expect_lte(max_type1(RR), alpha + 1e-10)
    }
  }
})

test_that("a finer nuisance parameter grid gives larger p-values", {
  # The 199 point grid contains every point of the 100 point grid
  for (tst in c('Z-pool', 'Boschloo')) {
    coarse <- attr(BinaryRR(7, 7, 0.05, tst, n.grid = 100), 'p.value')
    fine <- attr(BinaryRR(7, 7, 0.05, tst, n.grid = 199), 'p.value')
    expect_true(all(fine >= coarse - 1e-12), info = tst)
  }
})

test_that("the Boschloo test is at least as powerful as the Fisher test", {
  for (alt in c('greater', 'two.sided')) {
    fisher <- as_plain(BinaryRR(9, 9, 0.05, 'Fisher', alternative = alt))
    boschloo <- as_plain(BinaryRR(9, 9, 0.05, 'Boschloo', alternative = alt, n.grid = 200))
    expect_true(all(boschloo[fisher]), info = alt)
    expect_gte(sum(boschloo), sum(fisher))
  }
})

test_that("the Boschloo region contains the Fisher region under each two-sided rule", {
  for (ts in c('minlike', 'central', 'blaker')) {
    fisher <- as_plain(BinaryRR(12, 7, 0.05, 'Fisher', alternative = 'two.sided',
                                tsmethod = ts))
    boschloo <- as_plain(BinaryRR(12, 7, 0.05, 'Boschloo', alternative = 'two.sided',
                                  tsmethod = ts))
    expect_gt(sum(fisher), 0)
    expect_true(all(boschloo[fisher]), info = ts)
  }
})

test_that("the blaker Boschloo p-value can exceed the central Boschloo p-value", {
  # The two conventions order the outcomes differently, so the bound of the exact Fisher
  # p-value does not carry over: 2 of 4 against 0 of 3 gives about 0.305 and 0.219
  p <- vapply(c('blaker', 'central'), function(ts) {
    attr(BinaryRR(4, 3, 0.05, 'Boschloo', alternative = 'two.sided', tsmethod = ts),
         'p.value')[3, 1]
  }, numeric(1))
  expect_gt(p[['blaker']], p[['central']] + 0.05)
})

test_that("the blaker Boschloo test orders the outcomes by the blaker Fisher p-value", {
  N1 <- 5
  N2 <- 4
  p <- attr(BinaryRR(N1, N2, 0.05, 'Boschloo', alternative = 'two.sided',
                     tsmethod = 'blaker', n.grid = 25), 'p.value')
  ref <- unconditional_ref(fisher_blaker_ref(N1, N2), N1, N2, 25, decreasing = FALSE)
  expect_equal(matrix(as.vector(p), nrow = N1 + 1), ref, tolerance = 1e-12)
  # Sums of the p-values on the default grid and refined, and the numbers of outcomes
  # rejected at the level 0.05 by the Fisher test and by the refined Boschloo test, for 14
  # against 7 patients, from tools/reference/reference_values.py
  grid <- attr(BinaryRR(14, 7, 0.05, 'Boschloo', alternative = 'two.sided',
                        tsmethod = 'blaker'), 'p.value')
  expect_equal(sum(grid), 37.7030207902254, tolerance = 1e-10)
  ref <- BinaryRR(14, 7, 0.05, 'Boschloo', alternative = 'two.sided', tsmethod = 'blaker',
                  ref.pvalue = TRUE)
  expect_equal(sum(attr(ref, 'p.value')), 37.7046195635696, tolerance = 1e-10)
  expect_equal(sum(BinaryRR(14, 7, 0.05, 'Fisher', alternative = 'two.sided',
                            tsmethod = 'blaker')), 38)
  expect_equal(sum(ref), 44)
  seen <- ref_pvalue_calls(BinaryRR(14, 7, 0.05, 'Boschloo', alternative = 'two.sided',
                                    tsmethod = 'blaker'))
  expect_equal(seen$tsmethod, 'blaker')
  expect_equal(attr(ref, 'tsmethod'), 'blaker')
})

test_that("the Berger-Boos p-value matches a direct evaluation", {
  N1 <- 6
  N2 <- 5
  n.grid <- 25
  gam <- 0.001
  stat <- fisher_pvalue(N1, N2, 'greater', 'minlike', midp = FALSE)
  expect_equal(
    unconditional_pvalue(stat, N1, N2, n.grid, gam, decreasing = FALSE,
                         ref.pvalue = FALSE),
    berger_boos_ref(stat, N1, N2, n.grid, gam, decreasing = FALSE),
    tolerance = 1e-12
  )
  stat <- zstat(N1, N2)
  expect_equal(
    unconditional_pvalue(stat, N1, N2, n.grid, gam, decreasing = TRUE,
                         ref.pvalue = FALSE),
    berger_boos_ref(stat, N1, N2, n.grid, gam, decreasing = TRUE),
    tolerance = 1e-12
  )
})

test_that("the Berger-Boos procedure keeps the p-value above gamma and the level valid", {
  N1 <- 8
  N2 <- 8
  alpha <- 0.05
  gam <- 0.001
  for (tst in c('Z-pool', 'Boschloo')) {
    bb <- BinaryRR(N1, N2, alpha, tst, n.grid = 50, bb.gamma = gam)
    p.bb <- attr(bb, 'p.value')
    expect_true(all(p.bb >= gam - 1e-12), info = tst)
    expect_true(all(p.bb <= 1 + 1e-12), info = tst)
    expect_lte(max_type1(bb), alpha + 1e-10)
  }
})

test_that("Berger-Boos is ignored by the conditional tests, with a warning", {
  expect_warning(BinaryRR(5, 5, 0.05, 'Fisher', bb.gamma = 0.001), 'bb.gamma')
})

test_that("BinaryRR validates its arguments", {
  expect_error(BinaryRR(0, 5, 0.05, 'Chisq'), 'positive integers')
  expect_error(BinaryRR(5.5, 5, 0.05, 'Chisq'), 'positive integers')
  expect_error(BinaryRR(5, 5, 0, 'Chisq'), 'alpha')
  expect_error(BinaryRR(5, 5, 1, 'Chisq'), 'alpha')
  expect_error(BinaryRR(5, 5, 0.05, 'nonsense'), 'should be one of')
  expect_error(BinaryRR(5, 5, 0.05, 'Boschloo', n.grid = 1), 'n.grid')
  expect_error(BinaryRR(5, 5, 0.05, 'Boschloo', bb.gamma = -1), 'bb.gamma')
  expect_error(BinaryRR(5, 5, 0.05, 'Boschloo', bb.gamma = 0.05), 'bb.gamma')
  expect_error(BinaryRR(5, 5, 0.05, 'Chisq', alternative = 'lower'),
               'should be one of')
})

test_that("a smaller significance level gives a smaller rejection region", {
  for (tst in all_tests) {
    small <- as_plain(BinaryRR(8, 8, 0.01, tst, n.grid = 30))
    large <- as_plain(BinaryRR(8, 8, 0.05, tst, n.grid = 30))
    expect_true(all(large[small]), info = tst)
  }
})

test_that("the p-value attribute keeps the shape of the outcome grid", {
  # pmin() with a scalar first argument drops the dim attribute of its second argument
  for (tst in c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')) {
    for (alt in c('greater', 'two.sided')) {
      RR <- BinaryRR(4, 5, 0.05, tst, alternative = alt, n.grid = 20)
      expect_equal(dim(attr(RR, 'p.value')), c(5L, 6L),
                   info = sprintf('%s, %s', tst, alt))
    }
  }
})

test_that("the lower-tail region rejects only when group 1 falls behind", {
  N1 <- 9
  N2 <- 7
  rd <- outer(0:N1 / N1, 0:N2 / N2, '-')
  for (tst in all_tests) {
    RR <- BinaryRR(N1, N2, 0.05, tst, alternative = 'less', n.grid = 30)
    expect_equal(attr(RR, 'alternative'), 'less')
    expect_true(all(rd[as_plain(RR)] < 0), info = tst)
    expect_true(any(as_plain(RR)), info = tst)
  }
})

# The expected values below are certified by a branch and bound on the Bernstein
# coefficients of the tail probability, see tools/reference/reference_values.py

test_that("ref.pvalue removes the excess size caused by the grid maximum", {
  size_of <- function(RR) {
    f <- function(t) {
      vapply(t, function(u) power_from_rr(RR, dbinom(0:32, 32, u), dbinom(0:32, 32, u)),
             numeric(1))
    }
    x <- seq(0.005, 0.995, by = 0.005)
    refine_max(f, x, f(x))$y
  }
  rejected <- list('Z-pool' = c(342, 340), Boschloo = c(338, 336))
  size <- list('Z-pool' = c(0.0250058337172828, 0.023344357650963),
               Boschloo = c(0.0250047740156628, 0.0233381055506759))
  for (tst in c('Z-pool', 'Boschloo')) {
    grid <- BinaryRR(32, 32, 0.025, tst)
    ref <- BinaryRR(32, 32, 0.025, tst, ref.pvalue = TRUE)
    expect_true(attr(ref, 'ref.pvalue'))
    expect_equal(c(sum(grid), sum(ref)), rejected[[tst]], info = tst)
    expect_true(all(attr(ref, 'p.value') >= attr(grid, 'p.value')), info = tst)
    got <- c(size_of(as_plain(grid)), size_of(as_plain(ref)))
    expect_equal(got, size[[tst]], tolerance = 1e-9, info = tst)
    expect_gt(got[1], 0.025)
    expect_lte(got[2], 0.025)
  }
})

test_that("ref.pvalue resolves a maximum that lies between two uniform grid points", {
  grid <- attr(BinaryRR(150, 60, 0.025, 'Z-pool'), 'p.value')
  ref <- attr(BinaryRR(150, 60, 0.025, 'Z-pool', ref.pvalue = TRUE), 'p.value')
  expect_equal(ref[104, 33], 0.0304116803087155, tolerance = 1e-10)
  expect_true(all(ref >= grid))
})

test_that("ref.pvalue copes with grid points that coincide up to rounding", {
  # With n.grid = 101 the uniform grid contains 0.25, 0.5 and 0.75, which the arcsine grid
  # also contains up to rounding. The certified sums do not depend on the grid
  p <- attr(BinaryRR(6, 12, 0.025, 'Z-pool', n.grid = 101, ref.pvalue = TRUE), 'p.value')
  expect_equal(sum(p), 55.1616136367008, tolerance = 1e-10)
  p <- attr(BinaryRR(5, 89, 0.025, 'Boschloo', n.grid = 101, ref.pvalue = TRUE),
            'p.value')
  expect_equal(sum(p), 289.919122001934, tolerance = 1e-10)
})

test_that("ref.pvalue applies to the lower-tail and two-sided versions", {
  for (alt in c('less', 'two.sided')) {
    grid <- attr(BinaryRR(20, 14, 0.05, 'Z-pool', alternative = alt), 'p.value')
    ref <- attr(BinaryRR(20, 14, 0.05, 'Z-pool', alternative = alt, ref.pvalue = TRUE),
                'p.value')
    expect_true(all(ref >= grid), info = alt)
    expect_gt(max(ref - grid), 1e-6)
  }
  # The lower-tail version is the upper-tail version with the groups exchanged
  less <- attr(BinaryRR(20, 14, 0.05, 'Z-pool', alternative = 'less', ref.pvalue = TRUE),
               'p.value')
  greater <- attr(BinaryRR(14, 20, 0.05, 'Z-pool', ref.pvalue = TRUE), 'p.value')
  expect_equal(less, t(greater))
})

test_that("ref.pvalue is ignored by the conditional tests and validated", {
  for (tst in c('Chisq', 'Fisher', 'Fisher-midP')) {
    expect_identical(BinaryRR(12, 9, 0.05, tst, ref.pvalue = TRUE),
                     BinaryRR(12, 9, 0.05, tst), info = tst)
  }
  expect_error(BinaryRR(5, 5, 0.05, 'Z-pool', ref.pvalue = NA), 'ref.pvalue')
  expect_error(BinaryRR(5, 5, 0.05, 'Z-pool', ref.pvalue = 'yes'), 'ref.pvalue')
  expect_error(BinaryRR(5, 5, 0.05, 'Z-pool', ref.pvalue = c(TRUE, TRUE)), 'ref.pvalue')
})

test_that("BinaryRR passes ref.pvalue to the p-value computation", {
  seen <- ref_pvalue_calls(BinaryRR(8, 6, 0.05, 'Boschloo', ref.pvalue = TRUE))
  expect_equal(seen$ref.pvalue, TRUE)
  seen <- ref_pvalue_calls(BinaryRR(8, 6, 0.05, 'Boschloo'))
  expect_equal(seen$ref.pvalue, FALSE)
})

test_that("the Blackwelder test reproduces the example of Blackwelder (1982)", {
  # 18 of 30 responders with the experimental and 21 of 30 with the standard therapy give
  # a one-sided p-value of 0.21 for a difference of 0.2 (p. 350)
  p <- attr(BinaryRR(30, 30, 0.05, 'Blackwelder', margin = 0.2), 'p.value')
  expect_equal(round(p[19, 22], 2), 0.21)
  expect_equal(p[19, 22], stats::pnorm(0.816496580927726, lower.tail = FALSE),
               tolerance = 1e-10)
})

test_that("a margin is tested in the direction of the alternative", {
  for (tst in c('Blackwelder', 'Farrington-Manning')) {
    less <- attr(BinaryRR(14, 11, 0.05, tst, alternative = 'less', margin = 0.1),
                 'p.value')
    greater <- attr(BinaryRR(11, 14, 0.05, tst, margin = 0.1), 'p.value')
    expect_equal(less, t(greater), info = tst)
  }
  # Under 'less', H0: p1 - p2 >= margin is rejected for small p1 - p2. The p-value of 4 / 14
  # against 6 / 11 from the Blackwelder statistic written out directly
  pl <- attr(BinaryRR(14, 11, 0.05, 'Blackwelder', alternative = 'less', margin = 0.1),
             'p.value')
  se <- sqrt((4 / 14) * (10 / 14) / 14 + (6 / 11) * (5 / 11) / 11)
  expect_equal(pl[5, 7], pnorm((4 / 14 - 6 / 11 - 0.1) / se), tolerance = 1e-12)
  # The Farrington-Manning test without a margin is the chi-squared test
  expect_identical(as_plain(BinaryRR(15, 12, 0.025, 'Farrington-Manning')),
                   as_plain(BinaryRR(15, 12, 0.025, 'Chisq')))
  # A larger margin rejects more outcomes under the Blackwelder test, whose standard
  # error does not depend on the margin
  a <- as_plain(BinaryRR(20, 20, 0.025, 'Blackwelder', margin = 0.1))
  b <- as_plain(BinaryRR(20, 20, 0.025, 'Blackwelder', margin = 0.2))
  expect_true(all(b[a]))
  expect_gt(sum(b), sum(a))
})

test_that("a non-zero margin needs a non-inferiority test and a one-sided alternative", {
  expect_error(BinaryRR(10, 10, 0.05, 'Chisq', margin = 0.1), 'Blackwelder')
  expect_error(BinaryRR(10, 10, 0.05, 'Blackwelder', alternative = 'two.sided',
                        margin = 0.1), 'one-sided')
  expect_error(BinaryRR(10, 10, 0.05, 'Blackwelder', margin = 1), 'margin')
  expect_equal(attr(BinaryRR(10, 10, 0.05, 'Blackwelder', margin = 0.1), 'margin'), 0.1)
})

test_that("BinaryRR passes the margin to the p-value computation", {
  seen <- ref_pvalue_calls(BinaryRR(8, 6, 0.05, 'Farrington-Manning', margin = 0.15))
  expect_equal(seen$margin, 0.15)
  seen <- ref_pvalue_calls(BinaryRR(8, 6, 0.05, 'Farrington-Manning'))
  expect_equal(seen$margin, 0)
})

test_that("BinaryRR tests a margin on the scale of the risk ratio", {
  # Numbers of rejected outcomes from tools/reference/reference_values.py, which takes the
  # lower tail of the statistic for 'less' instead of exchanging the groups
  fm <- BinaryRR(40, 30, 0.025, 'Farrington-Manning', margin = 0.8, margin.scale = 'RR')
  bw <- BinaryRR(30, 40, 0.025, 'Blackwelder', alternative = 'less', margin = 1.25,
                 margin.scale = 'RR')
  expect_equal(c(sum(fm), sum(bw)), c(537, 543))
  expect_equal(attr(fm, 'margin'), 0.8)
  expect_identical(attr(fm, 'margin.scale'), 'RR')
  expect_identical(attr(BinaryRR(10, 10, 0.05, 'Chisq'), 'margin.scale'), 'RD')
  # The null hypothesis p1 / p2 >= 1.25 is p2 / p1 <= 0.8 with the groups exchanged
  less <- as_plain(BinaryRR(30, 40, 0.025, 'Farrington-Manning', alternative = 'less',
                            margin = 1.25, margin.scale = 'RR'))
  expect_identical(less, t(as_plain(fm)))
  # The p-value of 24 of 30 against 30 of 40 from the statistic written out directly
  p <- attr(bw, 'p.value')
  se <- sqrt((0.8 * 0.2) / 30 + 1.25^2 * (0.75 * 0.25) / 40)
  expect_equal(p[25, 31], pnorm((0.8 - 1.25 * 0.75) / se), tolerance = 1e-12)
})

test_that("a ratio margin needs a non-inferiority test, a one-sided alternative and R0 > 0", {
  expect_error(BinaryRR(10, 10, 0.05, 'Chisq', margin = 0.8, margin.scale = 'RR'),
               'Blackwelder')
  expect_error(BinaryRR(10, 10, 0.05, 'Blackwelder', alternative = 'two.sided',
                        margin = 0.8, margin.scale = 'RR'), 'one-sided')
  expect_error(BinaryRR(10, 10, 0.05, 'Blackwelder', margin.scale = 'RR'), 'positive')
  expect_error(BinaryRR(10, 10, 0.05, 'Blackwelder', margin = 0.8, margin.scale = 'OR'),
               'should be one of')
})

test_that("BinaryRR passes the scale of the margin to the p-value computation", {
  seen <- ref_pvalue_calls(BinaryRR(8, 6, 0.05, 'Farrington-Manning', margin = 0.8,
                                    margin.scale = 'RR'))
  expect_identical(seen$margin.scale, 'RR')
  expect_equal(seen$margin, 0.8)
  seen <- ref_pvalue_calls(BinaryRR(8, 6, 0.05, 'Farrington-Manning', margin = 0.15))
  expect_identical(seen$margin.scale, 'RD')
})
