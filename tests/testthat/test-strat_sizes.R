test_that("strat_sizes divides the nominal total between the strata and groups", {
  sz <- strat_sizes(c(10, 50), c(0.4, 0.6), c(1, 2))
  # Group 2: ceiling(c(0.4, 0.6) n / c(2, 3)); group 1: ceiling(r N2)
  expect_identical(sz$N2, matrix(c(2L, 10L, 2L, 10L), 2))
  expect_identical(sz$N1, matrix(c(2L, 10L, 4L, 20L), 2))
  expect_identical(dim(strat_sizes(7, c(0.5, 0.5), c(1, 1))$N1), c(1L, 2L))
})

test_that("strat_sizes gives totals of at least the nominal size that do not decrease", {
  n <- 1:200
  sz <- strat_sizes(n, c(0.3, 0.45, 0.25), c(1, 1.5, 2))
  tot <- rowSums(sz$N1) + rowSums(sz$N2)
  expect_true(all(tot >= n))
  expect_true(all(diff(sz$N1) >= 0) && all(diff(sz$N2) >= 0))
  expect_true(all(sz$N1 >= 1) && all(sz$N2 >= 1))
})
