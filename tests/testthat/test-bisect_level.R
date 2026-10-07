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

test_that("bisect_level passes the location of the last failing level as the probe", {
  calls <- NULL
  # The location reported for a level is the level itself, so a probe names the level at
  # which it was found
  assess <- function(des, probe) {
    calls <<- rbind(calls, c(level = des, probe = if (is.null(probe)) NA else probe))
    list(x = des, y = 1.6 * des, bound = NA_real_)
  }
  res <- bisect_level(identity, assess, NULL, 0.025, 1e-3)
  expect_lte(res$level, 0.025 / 1.6)
  expect_gt(res$level, 0.025 / 1.6 - 1e-3 * 0.025)
  n <- nrow(calls)
  fail <- 1.6 * calls[, 'level'] > 0.025
  expect_true(is.na(calls[1, 'probe']))
  expect_true(fail[1])
  expect_true(any(fail[-1]) && any(!fail[-1]))
  # Each later call receives the location of the last call that failed before it
  last <- cummax(ifelse(fail, seq_len(n), 0L))[-n]
  expect_equal(unname(calls[-1, 'probe']), unname(calls[last, 'level']))
})

test_that("bisect_level probes at the certified location after a certification fails", {
  calls <- NULL
  # The grid sees half of the rate and accepts every level, and the certification reports
  # the location 1 + level, so a probe above 1 comes from a certification
  assess <- function(des, probe) {
    calls <<- rbind(calls, c(level = des, probe = if (is.null(probe)) NA else probe))
    list(x = des, y = 0.8 * des, bound = NA_real_)
  }
  certify <- function(des, m) list(x = 1 + des, y = 1.6 * des, bound = 1.6 * des + 1e-13)
  res <- bisect_level(identity, assess, certify, 0.025, 1e-3)
  expect_lte(res$m$bound, 0.025)
  expect_gt(res$level, 0.025 / 1.6 - 1e-3 * 0.025)
  n <- nrow(calls)
  # After the nominal level the first bisection accepts every level, so its levels
  # increase, and the second bisection starts below the level found by the first
  k <- which(diff(calls[-1, 'level']) < 0)[1] + 2
  expect_false(is.na(k))
  expect_gt(k, 3)
  # The first bisection probes at the certified location of the nominal level
  expect_equal(unname(calls[2:(k - 1), 'probe']), rep(1.025, k - 2))
  # The second starts at the certified location of the level found by the first, and then
  # probes where the last level failed its certification
  lv <- calls[k:n, 'level']
  fail <- 1.6 * lv > 0.025
  expect_true(any(fail) && any(!fail))
  idx <- cummax(ifelse(fail, seq_along(lv), 0L))[-length(lv)]
  start <- 1 + calls[k - 1, 'level']
  want <- c(start, ifelse(idx > 0, 1 + lv[pmax(idx, 1L)], start))
  expect_equal(unname(calls[k:n, 'probe']), unname(want))
})
