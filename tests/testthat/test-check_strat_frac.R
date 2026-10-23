test_that("check_strat_frac rescales the fractions and recycles the ratio", {
  f <- check_strat_frac(c(20, 30), 2, 2)
  expect_equal(f$frac, c(0.4, 0.6))
  expect_equal(f$r, c(2, 2))
  expect_equal(check_strat_frac(c(1, 1, 2), c(1, 2, 3), 3)$r, c(1, 2, 3))
})

test_that("check_strat_frac rejects invalid fractions and ratios", {
  expect_error(check_strat_frac(c(0.5, 0.5), 1, 3), 'frac')
  expect_error(check_strat_frac(c(0.5, -0.5), 1, 2), 'frac')
  expect_error(check_strat_frac(c(0.5, 0.5), c(1, 2, 3), 2), 'r must')
  expect_error(check_strat_frac(c(0.5, 0.5), 0, 2), 'r must')
})
