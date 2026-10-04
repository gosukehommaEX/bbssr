# Reference values computed with an independent Python implementation (scipy)
test_that("ss_raw_n2 reproduces the normal approximation", {
  expect_equal(ss_raw_n2(0.42, 0.27, 1, 0.025, 0.8, 'greater', 'standard'),
               156.473696640835, tolerance = 1e-12)
  expect_equal(ss_raw_n2(0.25, 0.05, 1, 0.05, 0.8, 'two.sided', 'null.variance'),
               50.0366083064754, tolerance = 1e-12)
  expect_equal(ss_raw_n2(0.95, 0.75, 3, 0.05, 0.8, 'two.sided', 'null.variance'),
               23.5466392030473, tolerance = 1e-12)
  expect_equal(ss_raw_n2(0.6, 0.3, 2, 0.025, 0.9, 'greater', 'standard'),
               41.6637220158309, tolerance = 1e-12)
})

test_that("ss_raw_n2 gives the published sample sizes", {
  # Kieser (2020), Example 21.1: 2 x 157 patients
  expect_equal(ceiling(ss_raw_n2(0.42, 0.27, 1, 0.05, 0.8, 'two.sided', 'standard')), 157)
  # Friede and Kieser (2004), Table I with theta = 1: 102, 166 and 198 patients
  n <- vapply(c(0.05, 0.2, 0.4), function(p) {
    2 * ceiling(ss_raw_n2(p + 0.2, p, 1, 0.05, 0.8, 'two.sided', 'null.variance'))
  }, numeric(1))
  expect_equal(n, c(102, 166, 198))
})

test_that("the standard method equals the starting value of the exact search", {
  for (r in c(0.5, 1, 2)) {
    p1 <- 0.55
    p2 <- 0.3
    p <- (r * p1 + p2) / (1 + r)
    init <- (1 + 1 / r) / ((p1 - p2)^2) *
      (stats::qnorm(0.025) * sqrt(p * (1 - p)) +
         stats::qnorm(0.2) * sqrt((p1 * (1 - p1) / r + p2 * (1 - p2)) / (1 + 1 / r)))^2
    expect_equal(ss_raw_n2(p1, p2, r, 0.025, 0.8, 'greater', 'standard'), init,
                 tolerance = 1e-12, info = sprintf('r = %g', r))
  }
})

test_that("ss_raw_n2 is vectorized and symmetric in the direction of the effect", {
  a <- ss_raw_n2(c(0.4, 0.5), c(0.2, 0.2), 1, 0.025, 0.8, 'greater', 'standard')
  b <- ss_raw_n2(c(0.2, 0.2), c(0.4, 0.5), 1, 0.025, 0.8, 'less', 'standard')
  expect_length(a, 2)
  expect_equal(a, b)
})

test_that("probabilities outside the unit interval keep the assumed difference", {
  # Pooled rate 0.02 with a difference of 0.2: the recovered rates are 0.12 and -0.08
  got <- ss_raw_n2(0.12, -0.08, 1, 0.05, 0.8, 'two.sided', 'null.variance')
  z <- stats::qnorm(0.975) + stats::qnorm(0.8)
  expect_equal(got, z^2 * 2 * 0.02 * 0.98 / 0.2^2, tolerance = 1e-12)
  # Under the standard method the negative variance is truncated at zero
  got <- ss_raw_n2(0.12, -0.08, 1, 0.025, 0.8, 'greater', 'standard')
  want <- (stats::qnorm(0.975) * sqrt(2 * 0.02 * 0.98) +
             stats::qnorm(0.8) * sqrt(0.12 * 0.88))^2 / 0.2^2
  expect_equal(got, want, tolerance = 1e-12)
})
