test_that("strat_bern_eval evaluates a bivariate Bernstein polynomial", {
  # Constant coefficients give a constant polynomial
  expect_equal(strat_bern_eval(matrix(0.3, 4, 3), c(0, 0.4, 1), c(0.7, 0.2, 1)),
               rep(0.3, 3))
  # Degree one in each variable: the bilinear interpolation of the corners
  coef <- matrix(c(1, 2, 3, 5), 2)
  x <- 0.25
  y <- 0.6
  expect_equal(strat_bern_eval(coef, x, y),
               (1 - x) * (1 - y) * 1 + x * (1 - y) * 2 + (1 - x) * y * 3 + x * y * 5)
  expect_error(strat_bern_eval(coef, 1.2, 0.5), 'unit square')
  expect_error(strat_bern_eval(coef, c(0.2, 0.3), 0.5), 'same length')
})
