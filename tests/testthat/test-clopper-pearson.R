test_that("ClopperPearsonCI works correctly", {
  # Basic functionality test
  result <- ClopperPearsonCI(x = 7, n = 20, alpha = 0.05)

  expect_type(result, "double")
  expect_length(result, 2)
  expect_true(result[1] <= result[2])  # Lower bound <= Upper bound
  expect_true(result[1] >= 0 && result[2] <= 1)  # Within [0,1]
})

test_that("ClopperPearsonCI handles edge cases", {
  # Test with x = 0 (no successes)
  result_zero <- ClopperPearsonCI(x = 0, n = 20, alpha = 0.05)
  expect_equal(result_zero[1], 0)  # Lower bound should be exactly 0
  expect_true(result_zero[2] > 0)  # Upper bound should be > 0

  # Test with x = n (all successes)
  result_all <- ClopperPearsonCI(x = 20, n = 20, alpha = 0.05)
  expect_true(result_all[1] < 1)  # Lower bound should be < 1
  expect_equal(result_all[2], 1)  # Upper bound should be exactly 1

  # Test with n = 0 (edge case)
  result_n_zero <- ClopperPearsonCI(x = 0, n = 0, alpha = 0.05)
  expect_equal(result_n_zero, c(0, 1))  # Should return [0, 1]
})

test_that("ClopperPearsonCI validates inputs", {
  # Test invalid x values
  expect_error(ClopperPearsonCI(x = -1, n = 20, alpha = 0.05))
  expect_error(ClopperPearsonCI(x = 21, n = 20, alpha = 0.05))

  # Test invalid n values
  expect_error(ClopperPearsonCI(x = 5, n = -1, alpha = 0.05))

  # Test invalid alpha values
  expect_error(ClopperPearsonCI(x = 7, n = 20, alpha = 0))
  expect_error(ClopperPearsonCI(x = 7, n = 20, alpha = 1))
  expect_error(ClopperPearsonCI(x = 7, n = 20, alpha = -0.1))
  expect_error(ClopperPearsonCI(x = 7, n = 20, alpha = 1.1))
})

test_that("ClopperPearsonCI confidence levels work correctly", {
  # Test different confidence levels
  result_90 <- ClopperPearsonCI(x = 7, n = 20, alpha = 0.1)  # 90% CI
  result_95 <- ClopperPearsonCI(x = 7, n = 20, alpha = 0.05) # 95% CI
  result_99 <- ClopperPearsonCI(x = 7, n = 20, alpha = 0.01) # 99% CI

  # Higher confidence level should give wider intervals
  interval_90 <- diff(result_90)
  interval_95 <- diff(result_95)
  interval_99 <- diff(result_99)

  expect_true(interval_90 < interval_95)
  expect_true(interval_95 < interval_99)
})

test_that("ClopperPearsonCI is consistent with known values", {
  # Test with known exact values for small samples
  # For x=1, n=10, alpha=0.05, the CI should be approximately [0.0025, 0.445]
  result_known <- ClopperPearsonCI(x = 1, n = 10, alpha = 0.05)
  expect_true(abs(result_known[1] - 0.0025) < 0.001)
  expect_true(abs(result_known[2] - 0.445) < 0.01)

  # For x=5, n=10, alpha=0.05, the proportion is 0.5, CI should be roughly [0.187, 0.813]
  result_half <- ClopperPearsonCI(x = 5, n = 10, alpha = 0.05)
  expect_true(result_half[1] > 0.15 && result_half[1] < 0.25)
  expect_true(result_half[2] > 0.75 && result_half[2] < 0.85)
})

test_that("ClopperPearsonCI handles small samples", {
  # Test with very small sample sizes
  result_n1 <- ClopperPearsonCI(x = 0, n = 1, alpha = 0.05)
  expect_equal(result_n1[1], 0)
  expect_true(result_n1[2] < 1)

  result_n1_success <- ClopperPearsonCI(x = 1, n = 1, alpha = 0.05)
  expect_true(result_n1_success[1] > 0)
  expect_equal(result_n1_success[2], 1)
})

test_that("ClopperPearsonCI handles large samples", {
  # Test with larger sample sizes
  result_large <- ClopperPearsonCI(x = 50, n = 100, alpha = 0.05)

  expect_type(result_large, "double")
  expect_length(result_large, 2)
  expect_true(result_large[1] <= 0.5 && result_large[2] >= 0.5)  # Should contain true proportion

  # Interval should be narrower for larger samples
  result_small <- ClopperPearsonCI(x = 5, n = 10, alpha = 0.05)
  interval_large <- diff(result_large)
  interval_small <- diff(result_small)

  expect_true(interval_large < interval_small)
})

test_that("ClopperPearsonCI bounds are correct", {
  # Test that bounds are always within [0, 1]
  test_cases <- list(
    c(0, 5), c(1, 5), c(2, 5), c(3, 5), c(4, 5), c(5, 5),
    c(0, 20), c(10, 20), c(20, 20),
    c(0, 100), c(50, 100), c(100, 100)
  )

  for (case in test_cases) {
    result <- ClopperPearsonCI(x = case[1], n = case[2], alpha = 0.05)
    expect_true(result[1] >= 0, info = paste("x =", case[1], ", n =", case[2]))
    expect_true(result[2] <= 1, info = paste("x =", case[1], ", n =", case[2]))
    expect_true(result[1] <= result[2], info = paste("x =", case[1], ", n =", case[2]))
  }
})
