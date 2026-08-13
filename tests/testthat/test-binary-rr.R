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
    for (ts in c('minlike', 'central')) {
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

test_that("the Berger-Boos p-value matches a direct evaluation", {
  N1 <- 6
  N2 <- 5
  n.grid <- 25
  gam <- 0.001
  stat <- fisher_pvalue(N1, N2, 'greater', 'minlike', midp = FALSE)
  expect_equal(
    unconditional_pvalue(stat, N1, N2, n.grid, gam, decreasing = FALSE),
    berger_boos_ref(stat, N1, N2, n.grid, gam, decreasing = FALSE),
    tolerance = 1e-12
  )
  stat <- zstat(N1, N2)
  expect_equal(
    unconditional_pvalue(stat, N1, N2, n.grid, gam, decreasing = TRUE),
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
  expect_error(BinaryRR(5, 5, 0.05, 'Chisq', alternative = 'less'),
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
