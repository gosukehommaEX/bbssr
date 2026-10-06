# Design of test-bssr_reject.R: Delta.A = 0.3, N1 = N2 = 10, interim sizes 4 and 3
bernstein_design <- function() {
  map <- bssr_map(0.3, 10, 10, NULL, c(4, 3), 1, 0.025, 0.8, 'Chisq', FALSE, 'greater',
                  'minlike', 100L, 0, 'RD', 'standard', 'Chisq', 0.025, 'group', NULL,
                  NULL, FALSE, 0)
  st <- bssr_setup(map)
  rr.list <- lapply(seq_along(st$N1), function(k) {
    get_rr(st$N1[k], st$N2[k], 0.025, 'Chisq', 'greater', 'minlike', 100, 0, FALSE, 0)
  })
  list(st = st, rr.list = rr.list)
}

# Value at t of the polynomial with the Bernstein coefficients coef
bernstein_value <- function(t, coef) {
  sum(coef * dbinom(seq_along(coef) - 1, length(coef) - 1, t))
}

test_that("tie_bernstein gives the rejection probability over the unit interval", {
  d <- bernstein_design()
  coef <- tie_bernstein(tie_weights(d$st, d$rr.list), c(0, 1), c(0, 1))
  # Final sizes of different totals contribute polynomials of different degrees
  expect_gt(length(unique(d$st$N1 + d$st$N2)), 1)
  expect_length(coef, max(d$st$N1 + d$st$N2) + 1)
  expect_true(all(coef >= 0))
  t <- c(0, 0.13, 0.5, 0.92, 1)
  expect_equal(vapply(t, bernstein_value, numeric(1), coef = coef),
               bssr_reject(d$st, d$rr.list, t, t), tolerance = 1e-13)
})

test_that("tie_bernstein follows response probabilities that are linear on an interval", {
  d <- bernstein_design()
  w <- tie_weights(d$st, d$rr.list)
  t <- c(0, 0.3, 0.71, 1)
  # Common response probability from 0.15 to 0.85
  coef <- tie_bernstein(w, c(0.15, 0.85), c(0.15, 0.85))
  th <- 0.15 + 0.7 * t
  expect_equal(vapply(t, bernstein_value, numeric(1), coef = coef),
               bssr_reject(d$st, d$rr.list, th, th), tolerance = 1e-13)
  # Different response probabilities, as on the boundary of a non-inferiority hypothesis
  coef <- tie_bernstein(w, c(0.05, 0.6), c(0.25, 0.8))
  expect_true(all(coef >= 0))
  expect_equal(vapply(t, bernstein_value, numeric(1), coef = coef),
               bssr_reject(d$st, d$rr.list, 0.05 + 0.55 * t, 0.25 + 0.55 * t),
               tolerance = 1e-13)
})
