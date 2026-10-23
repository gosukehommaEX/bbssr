test_that("paired_n1 returns the number of interim pairs", {
  expect_identical(paired_n1(30, 0.5, NULL), 15L)
  expect_identical(paired_n1(31, 0.5, NULL), 16L)
  expect_identical(paired_n1(NULL, NULL, 12), 12L)
})

test_that("paired_n1 checks its arguments", {
  expect_error(paired_n1(30, 0.5, 10), 'not both')
  expect_error(paired_n1(30, NULL, NULL), 'either omega or n.interim')
  expect_error(paired_n1(30, 1.5, NULL), 'omega must be')
  expect_error(paired_n1(30, NULL, 2.5), 'n.interim must be')
  expect_error(paired_n1(NULL, 0.5, NULL), 'N must be')
})

test_that("paired_n1 is not raised by a rounding error of the product", {
  # 0.55 * 100 exceeds 55 by a rounding error
  expect_identical(paired_n1(100, 0.55, NULL), 55L)
})
