test_that("ceil_tol ignores excesses caused by floating-point error", {
  expect_equal(ceil_tol(0.1 * 3 * 10), 3)
  expect_equal(ceiling(0.1 * 3 * 10), 4)
  expect_equal(ceil_tol(c(2.000001, 2.5, -0.5)), c(3, 3, 0))
  expect_equal(ceil_tol(Inf), Inf)
})
