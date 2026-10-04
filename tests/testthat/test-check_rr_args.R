test_that("check_rr_args returns the arguments in computational form", {
  a <- check_rr_args(10, 8, 0.025, 'Bosch', 50, 0)
  expect_identical(a, list(N1 = 10L, N2 = 8L, n.grid = 50L, Test = 'Boschloo'))
})

test_that("check_rr_args rejects invalid arguments", {
  expect_error(check_rr_args(c(5, 6), 5, 0.05, 'Chisq', 100, 0), 'single value')
  expect_error(check_rr_args(0, 5, 0.05, 'Chisq', 100, 0), 'positive integers')
  expect_error(check_rr_args(5, 4.5, 0.05, 'Chisq', 100, 0), 'positive integers')
  expect_error(check_rr_args(5, 5, 0, 'Chisq', 100, 0), 'alpha')
  expect_error(check_rr_args(5, 5, 0.05, 'Chisq', 1, 0), 'n.grid')
  expect_error(check_rr_args(5, 5, 0.05, 'Z-pool', 100, 0.05), 'bb.gamma')
  expect_error(check_rr_args(5, 5, 0.05, 'Student', 100, 0), 'should be one of')
})

test_that("check_rr_args warns when bb.gamma is given for a conditional test", {
  expect_warning(check_rr_args(5, 5, 0.05, 'Fisher', 100, 0.001), 'bb.gamma')
  expect_no_warning(check_rr_args(5, 5, 0.05, 'Z-pool', 100, 0.001))
})
