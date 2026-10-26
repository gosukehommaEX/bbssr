test_that("PairedSampleSize reproduces the example of Connor (1987, p. 210)", {
  # theta1 = .80, theta2 = .70, one-sided level .05 and power .90: n = 254 for psi = .30
  # and n = 296 for psi = .35
  n <- function(p10, p01) {
    PairedSampleSize(p10, p01, 0.05, 0.9, 'McNemar', method = 'standard')$N
  }
  expect_equal(c(n(0.2, 0.1), n(0.225, 0.125)), c(254, 296))
})

test_that("PairedSampleSize reproduces the fixed designs of Schepers et al. (2026)", {
  # Smallest number of pairs at which the exact conditional McNemar test at the two-sided
  # level 0.05 has power 0.8: 72, 82 and 48 for p* = 0.875 and psi = 0.2, 0.175 and 0.3,
  # and 61, 92 and 112 for Delta = 0.15 and psi = 0.175, 0.25 and 0.3
  n <- function(psi, p10) {
    PairedSampleSize(p10, psi - p10, 0.05, 0.8, 'McNemar-exact', 'two.sided',
                     search = 'smallest')$N
  }
  expect_equal(c(n(0.2, 0.175), n(0.175, 0.153125), n(0.3, 0.2625)), c(72, 82, 48))
  expect_equal(c(n(0.175, 0.1625), n(0.25, 0.2), n(0.3, 0.225)), c(61, 92, 112))
})

test_that("PairedSampleSize agrees with the reference values", {
  std <- function(...) PairedSampleSize(..., method = 'standard')$N
  expect_equal(c(std(0.3, 0.1, 0.025, 0.8, 'McNemar'),
                 std(0.2, 0.1, 0.05, 0.9, 'McNemar', 'two.sided'),
                 std(0.1, 0.1, 0.025, 0.8, 'Paired-Score', margin = 0.1),
                 std(0.1, 0.25, 0.05, 0.8, 'McNemar', 'less')), c(77, 312, 171, 95))
  ex <- function(...) PairedSampleSize(...)$N
  expect_equal(c(ex(0.3, 0.1, 0.025, 0.8, 'McNemar-exact'),
                 ex(0.3, 0.1, 0.025, 0.8, 'McNemar'),
                 ex(0.1, 0.1, 0.025, 0.8, 'Paired-Score', margin = 0.1)), c(84, 72, 165))
  expect_equal(ex(0.3, 0.1, 0.025, 0.8, 'McNemar-midP', search = 'smallest'), 76)
})

test_that("PairedSampleSize of the continuity corrected test", {
  # The normal approximation does not depend on the test
  std <- function(Test) {
    PairedSampleSize(0.3, 0.1, 0.025, 0.8, Test, method = 'standard')$N
  }
  expect_equal(std('McNemar-cc'), std('McNemar'))
  expect_equal(c(PairedSampleSize(0.3, 0.1, 0.025, 0.8, 'McNemar-cc')$N,
                 PairedSampleSize(0.2, 0.1, 0.05, 0.9, 'McNemar-cc', 'two.sided')$N),
               c(84, 327))
})

test_that("PairedSampleSize returns the power at the number of pairs", {
  ss <- PairedSampleSize(0.3, 0.1, 0.025, 0.8, 'McNemar-exact')
  expect_s3_class(ss, 'bbssr_psamplesize')
  expect_named(ss, c('p10', 'p01', 'alpha', 'tar.power', 'Test', 'alternative', 'Power',
                     'N'))
  expect_equal(ss$Power, PairedPower(0.3, 0.1, ss$N, 0.025, 'McNemar-exact')$Power)
  expect_gte(ss$Power, 0.8)
  expect_equal(attr(ss, 'n.raw'), 76.0874844931802, tolerance = 1e-10)
})

test_that("PairedSampleSize checks the probabilities against the hypotheses", {
  expect_error(PairedSampleSize(0.1, 0.3, 0.025, 0.8, 'McNemar'), "exceed -margin")
  expect_error(PairedSampleSize(0.3, 0.1, 0.025, 0.8, 'McNemar', 'less'), 'fall below')
  expect_error(PairedSampleSize(0.2, 0.2, 0.05, 0.8, 'McNemar', 'two.sided'), 'must differ')
  expect_error(PairedSampleSize(0.4, 0.1, 0.05, 0.8, 'Paired-Score', 'equivalence',
                                margin = 0.2), 'strictly between')
  expect_error(PairedSampleSize(c(0.3, 0.4), 0.1, 0.025, 0.8, 'McNemar'), 'single value')
  expect_error(PairedSampleSize(0.3, 0.1, 0.025, 1, 'McNemar'), 'tar.power')
})
