# Tolerance of fpCompare, with which the package compares a p-value with a level
eps <- sqrt(.Machine$double.eps)

test_that("report_level moves the level found below the smallest p-value not rejected", {
  pv <- get_pvalue(39L, 39L, 'Chisq', 'greater', 'minlike', 100L, 0, FALSE, 0, 'RD')
  # Level found by the bisection of BinaryAlphaAdjBSSR for the fixed design with 39
  # patients per group, just below the smallest p-value not rejected plus eps
  level <- 0.020770048405393
  out <- report_level(list(pv), level)
  # Reference value from tools/reference/reference_values.py
  expect_equal(out, 0.02077, tolerance = 1e-12)
  rej <- pv %<<% level
  expect_gt(sum(rej), 0)
  # The level found rejects the smallest p-value not rejected when p <= level is rejected
  expect_false(identical(pv <= level, rej))
  # The reported level gives the same rejection region under the three rules
  expect_identical(pv %<<% out, rej)
  expect_identical(pv <= out, rej)
  expect_identical(pv < out, rej)
})

test_that("report_level uses more digits when six do not fit between the p-values", {
  p <- c(0.001, 0.01234567, 0.01234569, 0.05)
  level <- 0.01234569 + eps - 1e-10
  out <- report_level(list(matrix(p, 2, 2)), level)
  # Reference value from tools/reference/reference_values.py
  expect_equal(out, 0.012345689, tolerance = 1e-12)
  expect_identical(p <= out, p %<<% level)
})

test_that("report_level keeps the level when no other level gives the same region", {
  # The two smallest p-values are closer than eps
  p <- c(0.01, 0.01 + 1e-8, 0.03)
  level <- 0.01 + 1e-8 + eps - 1e-10
  expect_identical(report_level(list(p), level), level)
  # No p-value is rejected
  expect_identical(report_level(list(c(0.001, 0.02)), 0.0005), 0.0005)
})
