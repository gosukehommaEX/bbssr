test_that("paired_tie_weights sums the hypergeometric probabilities", {
  N.D1 <- c(8, 8, 10, 12, 12)
  w <- paired_tie_weights(N.D1, 4, c(8, 10, 12))
  expect_length(w, 3)
  expect_equal(lengths(w), c(9, 11, 13))
  expect_equal(w[[2]], stats::dhyper(2, 0:10, 10 - 0:10, 4))
  expect_equal(w[[1]], stats::dhyper(0, 0:8, 8 - 0:8, 4) + stats::dhyper(1, 0:8, 8 - 0:8, 4))
  # A single final number of pairs gives weight 1 for every number of discordant pairs
  expect_equal(paired_tie_weights(rep(9, 5), 4, 9)[[1]], rep(1, 10))
})
