test_that("bernstein_max finds the maximum of a Bernstein polynomial", {
  # 3 t^2 (1 - t), whose maximum 4 / 9 is attained at t = 2 / 3
  m <- bernstein_max(c(0, 0, 1, 0), 1e-12, 100000L)
  expect_lte(abs(m$value - 4 / 9), 1e-12)
  expect_gte(m$bound, 4 / 9)
  expect_lte(m$bound - m$value, 1e-12)
  expect_equal(m$x, 2 / 3, tolerance = 1e-5)
  expect_true(m$complete)
})

test_that("bernstein_max returns the end value when the maximum is at an end", {
  m <- bernstein_max(c(1, 0.5, 0.2), 1e-12, 100000L)
  expect_identical(c(m$x, m$value, m$bound), c(0, 1, 1))
  m <- bernstein_max(c(0.2, 0.5, 1), 1e-12, 100000L)
  expect_identical(c(m$x, m$value, m$bound), c(1, 1, 1))
  m <- bernstein_max(0.7, 1e-12, 100000L)
  expect_identical(c(m$value, m$bound), c(0.7, 0.7))
})

test_that("bernstein_max agrees with a fine grid and bounds the polynomial", {
  t <- seq(0, 1, length.out = 20001)
  basis <- outer(t, 0:30, function(t, k) dbinom(k, 30, t))
  set.seed(1)
  for (i in 1:5) {
    coef <- runif(31)
    m <- bernstein_max(coef, 1e-12, 100000L)
    v <- as.vector(basis %*% coef)
    expect_gte(m$value, max(v) - 1e-12)
    expect_gte(m$bound, m$value)
    expect_lte(m$bound - m$value, 1e-12)
    # The value is that of the polynomial at the reported location
    expect_equal(sum(dbinom(0:30, 30, m$x) * coef), m$value, tolerance = 1e-13)
  }
})

test_that("bernstein_max keeps a valid bound when the subdivision is cut short", {
  coef <- 1 + sin(seq(0, 3, length.out = 31))
  m <- bernstein_max(coef, 0, 2L)
  expect_false(m$complete)
  t <- seq(0, 1, length.out = 20001)
  v <- as.vector(outer(t, 0:30, function(t, k) dbinom(k, 30, t)) %*% coef)
  expect_gte(m$bound, max(v))
  expect_lte(m$value, max(v) + 1e-14)
})
