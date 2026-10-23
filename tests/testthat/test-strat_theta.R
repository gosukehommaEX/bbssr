test_that("strat_theta builds the points of the null hypothesis", {
  g <- strat_theta(c(0.5, 0.1, 0.9), 2)
  expect_identical(dim(g), c(9L, 2L))
  expect_identical(colnames(g), c('theta.S1', 'theta.S2'))
  expect_equal(g[1:4, 1], c(0.1, 0.5, 0.9, 0.1))
  expect_equal(g[1:4, 2], c(0.1, 0.1, 0.1, 0.5))
  d <- strat_theta(c(0.2, 0.4), 3)
  expect_equal(unname(d), cbind(c(0.2, 0.4), c(0.2, 0.4), c(0.2, 0.4)))
  m <- strat_theta(rbind(c(0.1, 0.2, 0.3)), 3)
  expect_equal(unname(m), rbind(c(0.1, 0.2, 0.3)))
})

test_that("strat_theta rejects invalid points", {
  expect_error(strat_theta(c(0.2, 1.2), 2), 'theta')
  expect_error(strat_theta(rbind(c(0.1, 0.2)), 3), 'one column per stratum')
})
