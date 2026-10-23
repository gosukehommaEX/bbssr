test_that("paired_raw_n reproduces Table 1 of Connor (1987)", {
  # Formula (3) with the quantiles 1.645 and 1.282 of the article
  za <- c(1.645, 1.645, 1.282, 1.282)
  zb <- c(1.282, 1.645, 1.282, 1.645)
  n <- c(mapply(function(a, b) {
    paired_raw_n(0.45, 0.05, stats::pnorm(-a), stats::pnorm(b), 'greater', 0)
  }, za, zb), mapply(function(a, b) {
    paired_raw_n(0.55, 0.05, stats::pnorm(-a), stats::pnorm(b), 'greater', 0)
  }, za, zb))
  expect_equal(round(n, 2), c(22.82, 28.15, 17.10, 21.76, 16.53, 20.20, 12.27, 15.46))
})

test_that("paired_raw_n reproduces Table 3 of Nam (1997)", {
  # n_t of formula (9) with z = 1.645 and 0.842, truncated to an integer. Two rows are
  # misprinted in the article: p10 = .05 in the row with p01 = .10 stands for .10, and
  # n_t = 188 in the last row stands for 118
  p01 <- rep(c(0.05, 0.10), c(8, 4))
  d0 <- rep(c(-0.03, -0.05, -0.05), each = 4)
  p10 <- c(0.05, 0.08, 0.10, 0.12, 0.05, 0.08, 0.10, 0.12, 0.10, 0.13, 0.15, 0.17)
  n <- mapply(function(p, q, d) {
    floor(paired_raw_n(p, q, stats::pnorm(-1.645), stats::pnorm(0.842), 'greater', -d))
  }, p10, p01, d0)
  expect_equal(n, c(717, 234, 151, 109, 273, 138, 100, 78, 506, 228, 158, 118))
})

test_that("paired_raw_n agrees with the reference values", {
  expect_equal(paired_raw_n(0.3, 0.1, 0.025, 0.8, 'greater', 0), 76.0874844931802,
               tolerance = 1e-10)
  expect_equal(paired_raw_n(0.1, 0.1, 0.025, 0.8, 'greater', 0.1), 170.132474774384,
               tolerance = 1e-10)
  expect_equal(paired_raw_n(0.12, 0.1, 0.05, 0.8, 'greater', -0.01), 13574.6280828215,
               tolerance = 1e-10)
  expect_equal(paired_raw_n(0.15, 0.15, 0.05, 0.8, 'equivalence', c(0.2, 0.2)),
               68.3571399132853, tolerance = 1e-10)
  expect_equal(paired_raw_n(0.2, 0.1, 0.05, 0.9, 'equivalence', c(0.2, 0.25)),
               115.160052104942, tolerance = 1e-10)
})

test_that("paired_raw_n uses the absolute difference for a two-sided alternative", {
  expect_equal(paired_raw_n(0.1, 0.2, 0.05, 0.9, 'two.sided', 0),
               paired_raw_n(0.2, 0.1, 0.025, 0.9, 'greater', 0))
})

test_that("paired_raw_n of a test of equivalence handles edge cases", {
  expect_true(is.nan(paired_raw_n(0.4, 0.1, 0.05, 0.8, 'equivalence', c(0.2, 0.2))))
  # With alpha above 0.5 the approximation exceeds the target without pairs
  expect_equal(paired_raw_n(0.15, 0.15, 0.9, 0.2, 'equivalence', c(0.2, 0.2)), 0)
  # Without discordant pairs the variance under the alternative is 0
  expect_true(is.finite(paired_raw_n(0, 0, 0.05, 0.8, 'equivalence', c(0.2, 0.2))))
})
