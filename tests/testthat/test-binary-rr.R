test_that("BinaryRR basic functionality", {
  # Basic test with small sample sizes and fastest test
  result <- BinaryRR(N1 = 5, N2 = 5, alpha = 0.05, Test = 'Chisq')

  expect_type(result, "logical")
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(6, 6))  # (N1+1) x (N2+1)
  expect_true(any(result == TRUE) || any(result == FALSE))  # Should contain logical values
})

test_that("BinaryRR matrix dimensions", {
  # Test matrix dimensions are correct
  N1 <- 3
  N2 <- 4
  result <- BinaryRR(N1 = N1, N2 = N2, alpha = 0.05, Test = 'Chisq')

  expect_equal(nrow(result), N1 + 1)
  expect_equal(ncol(result), N2 + 1)
  expect_true(is.logical(result))
})

test_that("BinaryRR different tests work", {
  # Test with different statistical tests using small samples
  tests <- c('Chisq', 'Fisher')

  for (test in tests) {
    result <- BinaryRR(N1 = 3, N2 = 3, alpha = 0.05, Test = test)

    expect_true(is.matrix(result))
    expect_true(is.logical(result))
    expect_equal(dim(result), c(4, 4))
  }
})

test_that("BinaryRR alpha levels", {
  # Test with different alpha levels
  result_05 <- BinaryRR(N1 = 4, N2 = 4, alpha = 0.05, Test = 'Chisq')
  result_01 <- BinaryRR(N1 = 4, N2 = 4, alpha = 0.01, Test = 'Chisq')

  # Smaller alpha should generally have fewer rejections
  rejections_05 <- sum(result_05)
  rejections_01 <- sum(result_01)

  expect_true(rejections_01 <= rejections_05)
})

test_that("BinaryRR edge cases", {
  # Test with minimal sample sizes
  result_min <- BinaryRR(N1 = 1, N2 = 1, alpha = 0.05, Test = 'Chisq')

  expect_equal(dim(result_min), c(2, 2))
  expect_true(is.logical(result_min))

  # Test with slightly larger but still small samples
  result_small <- BinaryRR(N1 = 2, N2 = 3, alpha = 0.1, Test = 'Chisq')

  expect_equal(dim(result_small), c(3, 4))
  expect_true(is.logical(result_small))
})

test_that("BinaryRR rejection patterns", {
  # Basic check that rejection region makes sense
  result <- BinaryRR(N1 = 4, N2 = 4, alpha = 0.05, Test = 'Chisq')

  # Should be logical matrix
  expect_true(is.logical(result))

  # Should have some structure (not all TRUE or all FALSE for reasonable alpha)
  expect_true(sum(result) > 0)  # At least some rejections
  expect_true(sum(result) < length(result))  # Not all rejections
})

test_that("BinaryRR parameter validation", {
  # Test that function handles basic parameter ranges
  # Using small values to avoid heavy computation

  # These should work without error
  expect_no_error(BinaryRR(N1 = 1, N2 = 1, alpha = 0.1, Test = 'Chisq'))
  expect_no_error(BinaryRR(N1 = 2, N2 = 2, alpha = 0.05, Test = 'Fisher'))
  expect_no_error(BinaryRR(N1 = 3, N2 = 3, alpha = 0.01, Test = 'Chisq'))
})
