test_that("strat_level_bracket brackets the level by the p-values of the design", {
  b <- strat_level_bracket(c(5L, 6L), c(6L, 4L), 0L, 0L, 0L, 0.05)
  expect_length(b, 2)
  expect_true(b[1] < 0.05 - 1.4901161193847656e-08 && b[2] >= 0.05 - 1.4901161193847656e-08)
  # Every level strictly between the two values, beyond the tolerance, gives the same
  # critical values
  a <- strat_crit(c(5L, 6L), c(6L, 4L), 0L, 0L, 0L, 0.05)
  mid <- (b[1] + 1.4901161193847656e-08 + b[2]) / 2
  m <- strat_crit(c(5L, 6L), c(6L, 4L), 0L, 0L, 0L, mid)
  expect_identical(m, a)
  # A slightly larger level rejects the p-value b[2]
  up <- strat_crit(c(5L, 6L), c(6L, 4L), 0L, 0L, 0L, b[2] + 1e-7)
  expect_false(identical(up, a))
  # Without a rejected p-value the lower value is -1
  expect_equal(strat_level_bracket(c(1L, 1L), c(1L, 1L), 4L, 0L, 0L, 0.01)[1], -1)
})
