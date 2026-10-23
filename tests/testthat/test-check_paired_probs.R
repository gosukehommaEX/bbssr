test_that("check_paired_probs accepts valid probabilities", {
  expect_null(check_paired_probs(c(0.2, 0.5), c(0.3, 0.5)))
  # A sum above 1 by rounding error is accepted
  expect_null(check_paired_probs(0.975, 0.025 + 1e-12))
})

test_that("check_paired_probs rejects invalid probabilities", {
  expect_error(check_paired_probs(0.2, c(0.1, 0.2)), 'same length')
  expect_error(check_paired_probs(-0.1, 0.2), 'lie in')
  expect_error(check_paired_probs(0.2, NA_real_), 'lie in')
  expect_error(check_paired_probs(0.6, 0.5), 'must not exceed 1')
})
