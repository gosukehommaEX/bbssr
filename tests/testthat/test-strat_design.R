test_that("strat_design derives the fractions and the interim sizes", {
  d <- strat_design(c(20, 30), c(20, 30), 1, NULL, 0.5, NULL, 0.025, 'CMH', 'greater',
                    'minlike', 'standard', 1e5)
  expect_equal(d$frac, c(0.4, 0.6))
  expect_identical(d$n1, c(10L, 15L))
  expect_identical(d$n2, c(10L, 15L))
  expect_identical(d$J, 2L)
  expect_equal(strat_design(c(20, 30), c(20, 30), 1, c(1, 1), 0.5, NULL, 0.025, 'CMH',
                            'greater', 'minlike', 'standard', 1e5)$frac, c(0.5, 0.5))
})

test_that("strat_design checks the settings", {
  expect_warning(strat_design(c(21, 30), c(20, 30), 1, NULL, 0.5, NULL, 0.025, 'CMH',
                              'greater', 'minlike', 'standard', 1e5), 'ceiling')
  expect_error(strat_design(c(20, 30, 10), c(20, 30, 10), 1, NULL, 0.5, NULL, 0.025, 'CMH',
                            'greater', 'minlike', 'exact', 1e5), 'standard')
  expect_error(strat_design(c(20, 30, 10), c(20, 30, 10), 1, NULL, 0.5, NULL, 0.025, 'CMH',
                            'greater', 'minlike', 'standard', 10), 'nsim')
})
