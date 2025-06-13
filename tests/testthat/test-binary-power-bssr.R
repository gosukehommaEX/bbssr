test_that("BinaryPowerBSSR basic functionality", {
  # Most basic test with minimal parameters
  result <- BinaryPowerBSSR(
    asmd.p1 = 0.9, asmd.p2 = 0.1,
    p = c(0.5), # Single pooled proportion only
    Delta.A = 0.8, Delta.T = 0.8,
    N1 = 2, N2 = 2, omega = 0.5, r = 1,
    alpha = 0.05, tar.power = 0.8,
    Test = 'Chisq',
    restricted = FALSE, weighted = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) == 1)
  expect_true("power.BSSR" %in% names(result))
  expect_true("power.TRAD" %in% names(result))
})

test_that("BinaryPowerBSSR output structure", {
  # Test output structure with minimal computation
  result <- BinaryPowerBSSR(
    asmd.p1 = 0.8, asmd.p2 = 0.2,
    p = c(0.4, 0.6), # Only 2 values
    Delta.A = 0.6, Delta.T = 0.6,
    N1 = 2, N2 = 2, omega = 0.5, r = 1,
    alpha = 0.05, tar.power = 0.8,
    Test = 'Chisq',
    restricted = FALSE, weighted = FALSE
  )

  expected_cols <- c("p1", "p2", "p", "power.BSSR", "power.TRAD")
  expect_true(all(expected_cols %in% names(result)))
  expect_equal(nrow(result), 2)  # Should match length of p vector
})

test_that("BinaryPowerBSSR power values are valid", {
  # Test that power values are within valid range
  # Use smaller Delta to ensure p1, p2 stay in valid range
  result <- BinaryPowerBSSR(
    asmd.p1 = 0.7, asmd.p2 = 0.3,
    p = c(0.4, 0.6),
    Delta.A = 0.4, Delta.T = 0.4,
    N1 = 2, N2 = 2, omega = 0.5, r = 1,
    alpha = 0.05, tar.power = 0.8,
    Test = 'Chisq',
    restricted = FALSE, weighted = FALSE
  )

  expect_true(all(result$power.BSSR >= 0 & result$power.BSSR <= 1))
  expect_true(all(result$power.TRAD >= 0 & result$power.TRAD <= 1))
})

test_that("BinaryPowerBSSR different design approaches", {
  # Test different BSSR approaches with minimal parameters
  designs <- list(
    list(restricted = TRUE, weighted = FALSE),
    list(restricted = FALSE, weighted = FALSE),
    list(restricted = FALSE, weighted = TRUE)
  )

  for (design in designs) {
    result <- BinaryPowerBSSR(
      asmd.p1 = 0.8, asmd.p2 = 0.2,
      p = c(0.5), # Single value only
      Delta.A = 0.6, Delta.T = 0.6,
      N1 = 2, N2 = 2, omega = 0.5, r = 1,
      alpha = 0.05, tar.power = 0.8,
      Test = 'Chisq',
      restricted = design$restricted,
      weighted = design$weighted
    )

    expect_s3_class(result, "data.frame")
    expect_true("power.BSSR" %in% names(result))
    expect_true(result$power.BSSR >= 0 && result$power.BSSR <= 1)
  }
})

test_that("BinaryPowerBSSR allocation ratios", {
  # Test different allocation ratios with safer parameters
  result_r1 <- BinaryPowerBSSR(
    asmd.p1 = 0.6, asmd.p2 = 0.4,
    p = c(0.5),
    Delta.A = 0.2, Delta.T = 0.2,
    N1 = 2, N2 = 2, omega = 0.5, r = 1,
    alpha = 0.05, tar.power = 0.8,
    Test = 'Chisq',
    restricted = FALSE, weighted = FALSE
  )

  expect_s3_class(result_r1, "data.frame")
  expect_true(result_r1$power.BSSR >= 0 && result_r1$power.BSSR <= 1)

  result_r2 <- BinaryPowerBSSR(
    asmd.p1 = 0.6, asmd.p2 = 0.4,
    p = c(0.5),
    Delta.A = 0.2, Delta.T = 0.2,
    N1 = 4, N2 = 2, omega = 0.5, r = 2,
    alpha = 0.05, tar.power = 0.8,
    Test = 'Chisq',
    restricted = FALSE, weighted = FALSE
  )

  expect_s3_class(result_r2, "data.frame")
  expect_true(result_r2$power.BSSR >= 0 && result_r2$power.BSSR <= 1)
})

test_that("BinaryPowerBSSR parameter preservation", {
  # Test that input parameters are preserved in output
  test_p <- c(0.4, 0.6)

  result <- BinaryPowerBSSR(
    asmd.p1 = 0.8, asmd.p2 = 0.2,
    p = test_p,
    Delta.A = 0.6, Delta.T = 0.6,
    N1 = 2, N2 = 2, omega = 0.5, r = 1,
    alpha = 0.05, tar.power = 0.8,
    Test = 'Chisq',
    restricted = FALSE, weighted = FALSE
  )

  expect_equal(result$p, test_p)
  expect_true(all(result$p1 == result$p + (1/(1+1)) * 0.6))  # Check p1 calculation
  expect_true(all(result$p2 == result$p - (1/(1+1)) * 0.6))  # Check p2 calculation
})

test_that("BinaryPowerBSSR minimal edge case", {
  # Test with absolute minimal parameters
  result <- BinaryPowerBSSR(
    asmd.p1 = 0.7, asmd.p2 = 0.3,
    p = 0.5, # Single value
    Delta.A = 0.4, Delta.T = 0.4,
    N1 = 1, N2 = 1, omega = 0.5, r = 1,
    alpha = 0.1, tar.power = 0.8,
    Test = 'Chisq',
    restricted = FALSE, weighted = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_type(result$power.BSSR, "double")
  expect_type(result$power.TRAD, "double")
})
