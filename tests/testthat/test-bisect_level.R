# Synthetic designs: a design is its level, the assessment on the grid sees the rate
# g(level), and the certification sees the rate cert(level), with an upper bound slightly
# above it
run_bisect <- function(g, cert, alpha = 0.025, tol = 1e-8) {
  assess <- function(des, probe) list(x = 0.5, y = g(des), bound = NA_real_)
  certify <- if (!is.null(cert)) {
    function(des, m) list(x = 0.3, y = cert(des), bound = cert(des) + 1e-13)
  }
  bisect_level(identity, assess, certify, alpha, tol)
}

test_that("bisect_level keeps the nominal level when it passes", {
  res <- run_bisect(function(l) l / 2, function(l) l / 2)
  expect_equal(res$level, 0.025)
  expect_identical(res$m, res$m0)
})

test_that("bisect_level finds the level at which the rate reaches alpha", {
  # The rate is 1.6 times the level, so the largest level that passes is alpha / 1.6
  res <- run_bisect(function(l) 1.6 * l, NULL)
  expect_lte(res$level, 0.025 / 1.6)
  expect_gt(res$level, 0.025 / 1.6 - 1e-8 * 0.025)
  expect_true(is.na(res$m$bound))
  res <- run_bisect(function(l) 1.6 * l, function(l) 1.6 * l)
  expect_equal(res$level, 0.025 / 1.6, tolerance = 1e-7)
  expect_lte(res$m$bound, 0.025)
})

test_that("bisect_level continues below a level whose certification fails", {
  # The grid sees half of the rate, so the bisection on the grid accepts every level and
  # the certification of the level it finds fails
  n.cert <- 0
  cert <- function(l) {
    n.cert <<- n.cert + 1
    1.6 * l
  }
  res <- run_bisect(function(l) 0.8 * l, cert)
  expect_equal(res$level, 0.025 / 1.6, tolerance = 1e-7)
  expect_lte(res$m$bound, 0.025)
  expect_equal(res$m0$y, 1.6 * 0.025)
  # Besides the nominal level and the level found on the grid, every level accepted
  # afterwards is certified
  expect_gt(n.cert, 2)
})

test_that("bisect_level returns the level 0 when no level passes", {
  res <- run_bisect(function(l) 1, NULL)
  expect_equal(res$level, 0)
  expect_equal(res$m$y, 0)
})
