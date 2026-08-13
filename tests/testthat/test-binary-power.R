all_tests <- c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')

test_that("BinaryPower returns a bbssr_power data frame", {
  res <- BinaryPower(p1 = 0.5, p2 = 0.2, N1 = 8, N2 = 8, alpha = 0.025, Test = 'Chisq')
  expect_s3_class(res, 'bbssr_power')
  expect_s3_class(res, 'data.frame')
  expect_equal(nrow(res), 1L)
  expect_named(res, c('p1', 'p2', 'N1', 'N2', 'alpha', 'Test', 'alternative', 'Power'))
  expect_type(res$Power, 'double')
  expect_true(res$Power >= 0 && res$Power <= 1)
})

test_that("BinaryPower is vectorised over the response probabilities", {
  res <- BinaryPower(p1 = c(0.4, 0.6, 0.8), p2 = rep(0.2, 3),
                     N1 = 10, N2 = 10, alpha = 0.025, Test = 'Chisq')
  expect_equal(nrow(res), 3L)
  expect_true(all(res$Power >= 0 & res$Power <= 1))
  expect_true(all(diff(res$Power) > 0))
})

test_that("BinaryPower agrees with a direct summation over the rejection region", {
  N1 <- 9
  N2 <- 6
  p1 <- 0.7
  p2 <- 0.25
  for (tst in all_tests) {
    RR <- BinaryRR(N1, N2, 0.025, tst, n.grid = 30)
    manual <- sum(outer(stats::dbinom(0:N1, N1, p1), stats::dbinom(0:N2, N2, p2)) *
                    as_plain(RR))
    got <- BinaryPower(p1, p2, N1, N2, 0.025, tst, n.grid = 30)$Power
    expect_equal(got, manual, tolerance = 1e-12, info = tst)
  }
})

test_that("BinaryPower rises with the sample size apart from the discreteness sawtooth", {
  # The exact power is not monotone in the sample size, because the attainable
  # significance level changes with each increment. At N1 = N2 = 8 the power of the
  # chi-squared test is 0.447 and at N1 = N2 = 10 it drops to 0.425, so only the overall
  # trend and a bound on the local decrease are asserted
  powers <- vapply(seq(6, 36, by = 6), function(n) {
    BinaryPower(0.6, 0.2, n, n, 0.025, 'Chisq')$Power
  }, numeric(1))
  expect_gt(powers[length(powers)], powers[1])
  expect_true(all(diff(powers) > -0.05))
})

test_that("BinaryPower is monotone in the effect size", {
  powers <- BinaryPower(seq(0.3, 0.8, by = 0.1), rep(0.2, 6), 15, 15, 0.025, 'Fisher')$Power
  expect_true(all(diff(powers) >= 0))
})

test_that("power under the null does not exceed the significance level", {
  alpha <- 0.05
  for (tst in c('Fisher', 'Z-pool', 'Boschloo')) {
    for (alt in c('greater', 'two.sided')) {
      p <- BinaryPower(0.4, 0.4, 10, 10, alpha, tst,
                       alternative = alt, n.grid = 100)$Power
      expect_lte(p, alpha + 1e-10)
    }
  }
})

test_that("the two-sided test is less powerful than the one-sided test", {
  one <- BinaryPower(0.7, 0.3, 15, 15, 0.025, 'Fisher')$Power
  two <- BinaryPower(0.7, 0.3, 15, 15, 0.025, 'Fisher', alternative = 'two.sided')$Power
  expect_lt(two, one)
})

test_that("BinaryPower validates its arguments", {
  expect_error(BinaryPower(c(0.4, 0.5), 0.2, 10, 10, 0.025, 'Chisq'), 'same length')
  expect_error(BinaryPower(-0.1, 0.2, 10, 10, 0.025, 'Chisq'), 'lie in')
  expect_error(BinaryPower(0.4, 1.2, 10, 10, 0.025, 'Chisq'), 'lie in')
})
