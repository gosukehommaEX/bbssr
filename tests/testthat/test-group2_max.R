test_that("group2_max agrees with the brute-force largest size of group 2", {
  # Expected values from a search over n2 with exact fractions for r; with r = 1.1,
  # 1.1 * 50 is slightly above 55 in floating point, and ceiling() alone gives 56
  got <- mapply(group2_max, c(1, 1.5, 2, 0.5, 1.1), c(50, 55, 31, 10, 105))
  expect_equal(got, c(25, 22, 10, 6, 50))
})

test_that("group2_max returns Inf without an upper bound and rejects an infinite bound", {
  expect_identical(group2_max(1, NULL), Inf)
  expect_error(group2_max(1, Inf), 'single finite number')
})
