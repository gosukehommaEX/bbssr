test_that("paired_cond_pvalue reproduces Table 6 of Fagerland et al. (2013)", {
  # Two-sided p-values of the exact and mid-p tests: Bentur et al. (N = 21, n12 = 1,
  # n21 = 7) and Cavo et al. (N = 161, n12 = 6, n21 = 16)
  expect_equal(round(paired_cond_pvalue(21, 'two.sided', FALSE)[2, 8], 4), 0.0703)
  expect_equal(round(paired_cond_pvalue(21, 'two.sided', TRUE)[2, 8], 4), 0.0391)
  expect_equal(round(paired_cond_pvalue(161, 'two.sided', FALSE)[7, 17], 4), 0.0525)
  expect_equal(round(paired_cond_pvalue(161, 'two.sided', TRUE)[7, 17], 4), 0.0347)
})

test_that("paired_cond_pvalue reproduces the conditional test of Berger and Sidik (2003)", {
  # Example of Section 4: x01 = 4, x10 = 11, n = 30, p-value of CB 0.0592
  expect_equal(round(paired_cond_pvalue(30, 'greater', FALSE)[12, 5], 4), 0.0592)
})

test_that("paired_cond_pvalue follows the binomial definitions", {
  P <- paired_cond_pvalue(12, 'greater', FALSE)
  M <- paired_cond_pvalue(12, 'greater', TRUE)
  T2 <- paired_cond_pvalue(12, 'two.sided', FALSE)
  T2m <- paired_cond_pvalue(12, 'two.sided', TRUE)
  expect_equal(P[6, 3], stats::pbinom(4, 7, 0.5, lower.tail = FALSE))
  expect_equal(M[6, 3], stats::pbinom(5, 7, 0.5, lower.tail = FALSE) +
                 stats::dbinom(5, 7, 0.5) / 2)
  expect_equal(T2[6, 3], 2 * stats::pbinom(2, 7, 0.5))
  expect_equal(T2m[6, 3], 2 * stats::pbinom(2, 7, 0.5) - stats::dbinom(5, 7, 0.5))
  # Tied counts: p-value 1 and mid-p value 1 - f / 2, formula (6) of Fagerland et al.
  expect_equal(T2[5, 5], 1)
  expect_equal(T2m[5, 5], 1 - stats::dbinom(4, 8, 0.5) / 2)
  # No discordant pairs
  expect_equal(c(P[1, 1], M[1, 1], T2[1, 1], T2m[1, 1]), c(1, 0.5, 1, 0.5))
  expect_true(all(is.na(P[outer(0:12, 0:12, '+') > 12])))
})

test_that("paired_cond_pvalue agrees with the reference values", {
  cells <- cbind(c(5, 6, 3, 0, 2, 4), c(1, 2, 0, 0, 3, 4)) + 1
  expect_equal(paired_cond_pvalue(12, 'greater', TRUE)[cells],
               c(0.0625, 0.08984375, 0.0625, 0.5, 0.65625, 0.5))
  expect_equal(paired_cond_pvalue(12, 'two.sided', TRUE)[cells],
               c(0.125, 0.1796875, 0.125, 0.5, 0.6875, 0.86328125))
})
