test_that("check_strat_probs returns a matrix with one column per stratum", {
  expect_identical(check_strat_probs(c(0.2, 0.4), 2, 'p1'), matrix(c(0.2, 0.4), 1))
  m <- rbind(c(0.2, 0.4, 0.6), c(0.1, 0.3, 0.5))
  expect_identical(check_strat_probs(m, 3, 'p1'), m)
})

test_that("check_strat_probs rejects invalid probabilities", {
  expect_error(check_strat_probs(c(0.2, 1.4), 2, 'p1'), 'p1 must contain')
  expect_error(check_strat_probs(c(0.2, 0.4, 0.5), 2, 'p2'), 'one element per stratum')
  expect_error(check_strat_probs(matrix(0.5, 2, 3), 2, 'p'), 'one column per stratum')
  expect_error(check_strat_probs(c(0.2, NA), 2, 'p'), 'p must contain')
})
