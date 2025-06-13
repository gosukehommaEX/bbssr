test_that("BinarySampleSize basic functionality", {
  # Most basic test with maximum effect size and fastest settings
  result <- BinarySampleSize(p1 = 1.0, p2 = 0.0, r = 1, alpha = 0.05, tar.power = 0.8, Test = 'Chisq')

  expect_s3_class(result, "data.frame")
  expect_true(result$N1 > 0)
  expect_true(result$N2 > 0)
  expect_equal(result$N, result$N1 + result$N2)
  expect_true(result$Power >= 0.8)
})

test_that("BinarySampleSize output structure", {
  # Test output structure without heavy computation
  result <- BinarySampleSize(p1 = 0.99, p2 = 0.01, r = 1, alpha = 0.05, tar.power = 0.8, Test = 'Chisq')

  expected_cols <- c("p1", "p2", "r", "alpha", "tar.power", "Test", "Power", "N1", "N2", "N")
  expect_true(all(expected_cols %in% names(result)))
  expect_equal(nrow(result), 1)
})

test_that("BinarySampleSize allocation ratio", {
  # Test allocation ratio with minimal computation
  result_r1 <- BinarySampleSize(p1 = 0.95, p2 = 0.05, r = 1, alpha = 0.05, tar.power = 0.8, Test = 'Chisq')
  expect_equal(result_r1$N1, result_r1$N2)  # Equal allocation
})

test_that("BinarySampleSize parameter preservation", {
  # Test that parameters are preserved in output
  result <- BinarySampleSize(p1 = 0.9, p2 = 0.1, r = 2, alpha = 0.05, tar.power = 0.8, Test = 'Chisq')

  expect_equal(result$p1, 0.9)
  expect_equal(result$p2, 0.1)
  expect_equal(result$r, 2)
  expect_equal(result$alpha, 0.05)
  expect_equal(result$tar.power, 0.8)
  expect_equal(result$Test, 'Chisq')
})
