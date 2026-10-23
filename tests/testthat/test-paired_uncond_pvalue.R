uncond_z <- function(N, alternative = 'greater', ...) {
  Z <- paired_score_z(N, 0)
  paired_uncond_pvalue(if (alternative == 'greater') Z else abs(Z), N, ...)
}

test_that("paired_uncond_pvalue reproduces the examples of Berger and Sidik (2003)", {
  # Section 4: x01 = 4, x10 = 11, n = 30, p-values of Z 0.0494 and of Z_C 0.0382 with
  # beta = 0.0005
  expect_equal(round(uncond_z(30, n.grid = 100, bb.gamma = 0, ref.pvalue = TRUE)[12, 5],
                     4), 0.0494)
  expect_equal(round(uncond_z(30, n.grid = 100, bb.gamma = 5e-4, ref.pvalue = TRUE)[12, 5],
                     4), 0.0382)
  # Second illustration: x01 = 2, x10 = 10, n = 20, p-values 0.01188 and 0.01238
  expect_equal(round(uncond_z(20, n.grid = 100, bb.gamma = 0, ref.pvalue = TRUE)[11, 3],
                     5), 0.01188)
  expect_equal(round(uncond_z(20, n.grid = 100, bb.gamma = 5e-4, ref.pvalue = TRUE)[11, 3],
                     5), 0.01238)
})

test_that("paired_uncond_pvalue reproduces Table 6 of Fagerland et al. (2013)", {
  # Bentur et al.: N = 21, n12 = 1, n21 = 7; the two-sided exact unconditional p-value of
  # the article is twice the one-sided p-value, 0.0353
  P <- t(uncond_z(21, n.grid = 100, bb.gamma = 0, ref.pvalue = TRUE))
  expect_equal(round(2 * P[2, 8], 4), 0.0353)
})

test_that("paired_uncond_pvalue agrees with the reference values", {
  cells <- cbind(c(5, 6, 3, 0, 2, 4), c(1, 2, 0, 0, 3, 4)) + 1
  expect_equal(uncond_z(12, n.grid = 100, bb.gamma = 0, ref.pvalue = TRUE)[cells],
               c(0.072998046875, 0.113135472540419, 0.072998046875, 1, 1, 1),
               tolerance = 1e-9)
  expect_equal(uncond_z(12, 'two.sided', n.grid = 100, bb.gamma = 0,
                        ref.pvalue = TRUE)[cells],
               c(0.14599609375, 0.226270945080838, 0.14599609375, 1, 0.7744140625, 1),
               tolerance = 1e-9)
  expect_equal(uncond_z(12, n.grid = 100, bb.gamma = 0.001, ref.pvalue = TRUE)[cells],
               c(0.0682777946325789, 0.114135472540419, 0.0539247380865849, 1,
                 0.767453288352893, 0.635274702815146), tolerance = 1e-9)
})

test_that("paired_uncond_pvalue gives outcomes with the same statistic one p-value", {
  P <- uncond_z(15, n.grid = 50, bb.gamma = 0, ref.pvalue = FALSE)
  Z <- paired_score_z(15, 0)
  ok <- !is.na(Z)
  # Ties exist, for example (2, 0) and (6, 2) both have Z = sqrt(2)
  expect_equal(Z[3, 1], Z[7, 3])
  expect_true(any(table(round(Z[ok], 9)) > 1))
  for (z in unique(round(Z[ok], 9))) {
    idx <- ok & abs(Z - z) < 1e-8
    expect_length(unique(round(P[idx], 13)), 1)
  }
})

test_that("refinement never lowers the p-value on the grid", {
  P0 <- uncond_z(18, n.grid = 30, bb.gamma = 0, ref.pvalue = FALSE)
  P1 <- uncond_z(18, n.grid = 30, bb.gamma = 0, ref.pvalue = TRUE)
  ok <- !is.na(P0)
  expect_true(all(P1[ok] >= P0[ok] - 1e-15))
})
