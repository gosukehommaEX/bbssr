test_that("strat_fixed_power uses the stored critical values with two strata", {
  p1 <- rbind(c(0.6, 0.45))
  p2 <- rbind(c(0.35, 0.3))
  expect_equal(strat_fixed_power(c(8, 6), c(7, 9), p1, p2, 0.05, 'Cochran-cc', 'greater',
                                 'minlike'),
               strat_power(c(8L, 6L), c(7L, 9L), p1, p2, 3L, 0L, 0L, 0.05),
               tolerance = 1e-13)
  # Three strata go through strat_power()
  p1 <- rbind(c(0.7, 0.3, 0.55))
  p2 <- rbind(c(0.4, 0.25, 0.5))
  expect_identical(strat_fixed_power(c(4, 5, 3), c(5, 3, 4), p1, p2, 0.1, 'CMH', 'less',
                                     'minlike'),
                   strat_power(c(4L, 5L, 3L), c(5L, 3L, 4L), p1, p2, 0L, 1L, 0L, 0.1))
})
