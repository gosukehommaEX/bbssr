test_that("ss_raw_parts gives the ingredients of ss_raw_n2", {
  rebuild <- function(p1, p2, r, alpha, tp, alt, method, margin, scale) {
    s <- ss_raw_parts(p1, p2, r, alt, method, margin, scale)
    a <- if (alt == 'two.sided') alpha / 2 else alpha
    (qnorm(1 - a) * sqrt(s$v0) + qnorm(tp) * sqrt(s$v1))^2 / s$d^2
  }
  cases <- list(
    list(0.6, 0.4, 1, 0.025, 0.8, 'greater', 'standard', 0, 'RD'),
    list(0.3, 0.45, 2, 0.05, 0.9, 'two.sided', 'null.variance', 0, 'RD'),
    list(0.7, 0.7, 3, 0.025, 0.8, 'greater', 'standard', 0.1, 'RD'),
    list(0.5, 0.55, 0.5, 0.025, 0.8, 'less', 'alternative.variance', 0.15, 'RD'),
    list(0.7, 0.7, 1, 0.025, 0.8, 'greater', 'standard', 0.8, 'RR'),
    list(0.5, 0.5, 1, 0.025, 0.8, 'less', 'standard', 1.25, 'RR')
  )
  for (cs in cases) {
    expect_equal(do.call(rebuild, cs), do.call(ss_raw_n2, cs), tolerance = 1e-12)
  }
})

test_that("ss_raw_parts gives the variances of the pooled Z test", {
  # p = 0.5, so v0 = 0.25 * 2 and v1 = 0.6 * 0.4 + 0.4 * 0.6
  s <- ss_raw_parts(0.6, 0.4, 1, 'greater', 'standard', 0, 'RD')
  expect_equal(c(s$d, s$v0, s$v1), c(0.2, 0.5, 0.48))
  s <- ss_raw_parts(0.6, 0.4, 1, 'greater', 'null.variance', 0, 'RD')
  expect_equal(c(s$v0, s$v1), c(0.5, 0.5))
})

test_that("ss_raw_parts measures the distance from the null boundary", {
  expect_equal(ss_raw_parts(0.7, 0.6, 1, 'greater', 'standard', 0.1, 'RD')$d, 0.2)
  expect_equal(ss_raw_parts(0.7, 0.6, 1, 'less', 'standard', 0.1, 'RD')$d, 0)
  expect_equal(ss_raw_parts(0.5, 0.5, 1, 'less', 'standard', 1.25, 'RR')$d, 0.125)
  # Under the alternative variance both variances are those of the true probabilities
  s <- ss_raw_parts(0.3, 0.4, 2, 'greater', 'alternative.variance', 0.1, 'RD')
  expect_equal(c(s$v0, s$v1), rep(0.3 * 0.7 / 2 + 0.4 * 0.6, 2))
})
