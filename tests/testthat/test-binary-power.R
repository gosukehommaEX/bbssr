test_that("BinaryPower works correctly", {
  # Basic functionality test
  result <- BinaryPower(p1 = 0.5, p2 = 0.3, N1 = 30, N2 = 30, alpha = 0.025, Test = 'Fisher')

  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result >= 0 && result <= 1)
})

test_that("BinaryPower handles edge cases", {
  # Test with equal probabilities (should give Type I error rate)
  result_null <- BinaryPower(p1 = 0.3, p2 = 0.3, N1 = 30, N2 = 30, alpha = 0.025, Test = 'Fisher')
  expect_true(result_null <= 0.03)  # Allow some tolerance for exact tests

  # Test with extreme probabilities
  result_extreme <- BinaryPower(p1 = 0.8, p2 = 0.2, N1 = 10, N2 = 10, alpha = 0.025, Test = 'Fisher')
  expect_true(result_extreme >= 0 && result_extreme <= 1)
})

test_that("BinaryPower works with vectors", {
  # Test with vector inputs
  p1_vec <- c(0.4, 0.5, 0.6)
  p2_vec <- c(0.2, 0.2, 0.2)

  results <- BinaryPower(p1 = p1_vec, p2 = p2_vec, N1 = 30, N2 = 30, alpha = 0.025, Test = 'Fisher')

  expect_length(results, 3)
  expect_true(all(results >= 0 & results <= 1))
  expect_true(results[3] > results[2])  # Higher effect size should have higher power
  expect_true(results[2] > results[1])
})

test_that("BinaryPower validates inputs", {
  # Test invalid probability values - expect warnings but not errors
  expect_warning(BinaryPower(p1 = -0.1, p2 = 0.3, N1 = 30, N2 = 30, alpha = 0.025, Test = 'Fisher'))
  expect_warning(BinaryPower(p1 = 1.1, p2 = 0.3, N1 = 30, N2 = 30, alpha = 0.025, Test = 'Fisher'))

  # Test small sample sizes - these should work
  result_small <- BinaryPower(p1 = 0.5, p2 = 0.3, N1 = 1, N2 = 30, alpha = 0.025, Test = 'Fisher')
  expect_type(result_small, "double")

  # Test edge alpha values
  result_alpha_small <- BinaryPower(p1 = 0.5, p2 = 0.3, N1 = 30, N2 = 30, alpha = 0.001, Test = 'Fisher')
  expect_type(result_alpha_small, "double")

  result_alpha_large <- BinaryPower(p1 = 0.5, p2 = 0.3, N1 = 30, N2 = 30, alpha = 0.5, Test = 'Fisher')
  expect_type(result_alpha_large, "double")

  # Test invalid test name - this should cause an error
  expect_error(BinaryPower(p1 = 0.5, p2 = 0.3, N1 = 30, N2 = 30, alpha = 0.025, Test = 'InvalidTest'))
})

test_that("All statistical tests work", {
  tests <- c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')

  for (test in tests) {
    result <- BinaryPower(p1 = 0.5, p2 = 0.3, N1 = 30, N2 = 30, alpha = 0.025, Test = test)
    expect_type(result, "double")
    expect_true(result >= 0 && result <= 1)
  }
})

test_that("BinaryPower is monotonic in effect size", {
  # Power should increase with effect size
  p2 <- 0.2
  p1_values <- seq(0.3, 0.8, by = 0.1)

  powers <- sapply(p1_values, function(p1) {
    BinaryPower(p1 = p1, p2 = p2, N1 = 30, N2 = 30, alpha = 0.025, Test = 'Fisher')
  })

  # Check monotonicity
  expect_true(all(diff(powers) >= 0))
})

test_that("BinaryPower is monotonic in sample size", {
  # Power should increase with sample size
  sample_sizes <- seq(10, 50, by = 10)

  powers <- sapply(sample_sizes, function(n) {
    BinaryPower(p1 = 0.5, p2 = 0.3, N1 = n, N2 = n, alpha = 0.025, Test = 'Fisher')
  })

  # Check monotonicity
  expect_true(all(diff(powers) >= 0))
})
