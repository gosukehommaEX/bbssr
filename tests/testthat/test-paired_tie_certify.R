test_that("paired_tie_certify bounds the largest rejection probability", {
  rr <- get_paired_rr(30, 0.05, 'McNemar', 'greater', 100, 0, FALSE, 0)
  cm <- paired_tie_certify(list(rr), 30, list(rep(1, 31)), c(0, 1), 'greater', 0)
  # Berger and Sidik (2003, Table 4): exact size 0.0558 of M with n = 30
  expect_equal(round(cm$y, 4), 0.0558)
  expect_equal(cm$y, 0.0558117064555376, tolerance = 1e-9)
  expect_gte(cm$bound, cm$y)
  expect_lte(cm$bound - cm$y, 1e-12 + 1e-15)
  b <- paired_boundary(cm$x, 'greater', 0)
  expect_equal(paired_power(rr, b$p10, b$p01), cm$y, tolerance = 1e-12)
  # A grid of values never exceeds the certified maximum
  psi <- seq(0, 1, by = 0.001)
  g <- paired_boundary(psi, 'greater', 0)
  expect_true(all(paired_power(rr, g$p10, g$p01) <= cm$bound + 1e-15))
})

test_that("paired_tie_certify works on a subinterval of a margin boundary", {
  rr <- get_paired_rr(25, 0.05, 'Paired-Score', 'greater', 100, 0, FALSE, 0.1)
  cm <- paired_tie_certify(list(rr), 25, list(rep(1, 26)), c(0.3, 0.8), 'greater', 0.1)
  expect_gte(cm$x, 0.3)
  expect_lte(cm$x, 0.8)
  psi <- seq(0.3, 0.8, by = 0.001)
  g <- paired_boundary(psi, 'greater', 0.1)
  vals <- paired_power(rr, g$p10, g$p01)
  expect_true(all(vals <= cm$bound + 1e-15))
  expect_gte(cm$y, max(vals) - 1e-15)
})
