test_that("strat_bssr_power without an interim analysis is the power of the fixed design", {
  crit <- get_strat_crit(c(7, 5), c(6, 8), 0.05, 'Strat-Fisher', 'two.sided', 'minlike')
  p1 <- rbind(c(0.7, 0.4), c(0.3, 0.3))
  p2 <- rbind(c(0.2, 0.3), c(0.3, 0.3))
  expect_equal(strat_bssr_power(c(0L, 0L), c(0L, 0L), list(crit), 0L, p1, p2),
               strat_power(c(7L, 5L), c(6L, 8L), p1, p2, 4L, 2L, 0L, 0.05),
               tolerance = 1e-13)
})

test_that("strat_bssr_power with a single final design is the power of that design", {
  # Whatever the interim outcome, the trial continues to the same final sizes
  crit <- get_strat_crit(c(7, 5), c(6, 8), 0.05, 'CMH', 'greater', 'minlike')
  id <- rep(0L, (3 + 2 + 1) * (2 + 4 + 1))
  p1 <- rbind(c(0.6, 0.5))
  p2 <- rbind(c(0.3, 0.2))
  expect_equal(strat_bssr_power(c(3L, 2L), c(2L, 4L), list(crit), id, p1, p2),
               strat_power(c(7L, 5L), c(6L, 8L), p1, p2, 0L, 0L, 0L, 0.05),
               tolerance = 1e-13)
})

test_that("strat_bssr_power checks its arguments", {
  crit <- get_strat_crit(c(7, 5), c(6, 8), 0.05, 'CMH', 'greater', 'minlike')
  p <- rbind(c(0.5, 0.5))
  expect_error(strat_bssr_power(c(3L, 2L), c(2L, 4L), list(crit), 0L, p, p),
               'one element per interim cell')
  expect_error(strat_bssr_power(c(3L, 2L), c(2L, 4L), list(crit), rep(1L, 42), p, p),
               'zero based')
  expect_error(strat_bssr_power(c(8L, 2L), c(2L, 4L), list(crit), rep(0L, 77), p, p),
               'below the interim')
  expect_error(strat_bssr_power(c(0L, 0L), c(0L, 0L), list(crit), 0L, p + 1, p),
               'unit interval|\\[0, 1\\]')
})
