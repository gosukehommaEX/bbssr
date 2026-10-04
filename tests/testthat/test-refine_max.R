test_that("refine_max finds a maximum between grid points", {
  f <- function(x) -(x - 0.3333)^2
  x <- seq(0, 1, by = 0.1)
  m <- refine_max(f, x, f(x))
  expect_equal(m$x, 0.3333, tolerance = 1e-6)
  expect_equal(m$y, 0, tolerance = 1e-10)
})

test_that("refine_max examines the largest local maxima", {
  # Two peaks: the grid favours the lower one, the refinement finds the higher one
  f <- function(x) pmax(1 - 200 * (x - 0.25)^2, 1.02 - 2000 * (x - 0.705)^2)
  x <- seq(0, 1, by = 0.05)
  expect_equal(x[which.max(f(x))], 0.25)
  m <- refine_max(f, x, f(x))
  expect_equal(m$x, 0.705, tolerance = 1e-6)
  expect_equal(m$y, 1.02, tolerance = 1e-10)
  # Without refinement the grid maximum is returned
  expect_lt(max(f(x)), 1.02)
})
