test_that("strat_bounds finds the range of the nominal size", {
  # Two equal strata with r = 1: every group has ceiling(n / 4) patients
  expect_equal(strat_bounds(c(0.5, 0.5), c(1, 1), c(5, 3), c(5, 3), NULL, NULL), c(17, Inf))
  expect_equal(strat_bounds(c(0.5, 0.5), c(1, 1), c(5, 3), c(5, 3), 30, 50), c(29, 48))
  expect_equal(strat_bounds(c(0.5, 0.5), c(1, 1), c(1, 1), c(1, 1), NULL, NULL), c(1, Inf))
})

test_that("strat_bounds rejects bounds that leave no design", {
  expect_error(strat_bounds(c(0.5, 0.5), c(1, 1), c(5, 3), c(5, 3), NULL, 16), 'N.max')
  expect_error(strat_bounds(c(0.5, 0.5), c(1, 1), c(1, 1), c(1, 1), c(1, 2), NULL),
               'N.min')
})
