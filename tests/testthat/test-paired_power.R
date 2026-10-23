test_that("paired_power sums the multinomial probabilities of the region", {
  set.seed(1)
  R <- random_region(9)
  for (p in list(c(0.2, 0.1), c(0.35, 0.35), c(0, 0.3), c(0.5, 0.5), c(0, 0))) {
    expect_equal(paired_power(R, p[1], p[2]), sum(trinomial_ref(9, p[1], p[2]) * R),
                 tolerance = 1e-13, info = paste(p, collapse = ', '))
  }
  # Vectorized over the probabilities
  expect_equal(paired_power(R, c(0.2, 0.35), c(0.1, 0.35)),
               c(sum(trinomial_ref(9, 0.2, 0.1) * R), sum(trinomial_ref(9, 0.35, 0.35) * R)),
               tolerance = 1e-13)
})

test_that("paired_power gives 1 for the region of all outcomes", {
  R <- outer(0:7, 0:7, '+') <= 7
  expect_equal(paired_power(R, c(0.1, 0.3, 0.6), c(0.2, 0.3, 0.4)), rep(1, 3),
               tolerance = 1e-13)
})

test_that("paired_power checks its arguments", {
  expect_error(paired_power(matrix(TRUE, 3, 4), 0.1, 0.1), 'square')
  expect_error(paired_power(matrix(TRUE, 3, 3), c(0.1, 0.2), 0.1), 'same length')
})
