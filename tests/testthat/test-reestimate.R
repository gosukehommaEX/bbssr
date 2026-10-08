test_that("reestimate evaluates each distinct pair once and reproduces BinarySampleSize", {
  hat.p1 <- c(0.5, 0.6, 0.5)
  hat.p2 <- c(0.2, 0.3, 0.2)
  re <- reestimate(hat.p1, hat.p2, 1, 0.025, 0.8, 'Fisher', 'greater', 'minlike', 100L,
                   0, 'exact', 'group', NULL, NULL, FALSE, 0, margin.scale = 'RD')
  expect_named(re, c('n.raw', 'N1.re', 'N2.re', 'N2.limit'))
  expect_true(all(is.na(re$N2.limit)))
  expect_true(all(is.na(re$n.raw)))
  for (i in 1:2) {
    ss <- BinarySampleSize(hat.p1[i], hat.p2[i], 1, 0.025, 0.8, 'Fisher')
    expect_equal(c(re$N1.re[i], re$N2.re[i]), c(ss$N1, ss$N2))
  }
  expect_equal(re[3, ], re[1, ], ignore_attr = TRUE)
})

test_that("the normal methods return the unrounded total", {
  re <- reestimate(0.4, 0.2, 2, 0.025, 0.8, 'Chisq', 'greater', 'minlike', 100L, 0,
                   'standard', 'group', NULL, NULL, FALSE, 0, margin.scale = 'RD')
  n2 <- ss_raw_n2(0.4, 0.2, 2, 0.025, 0.8, 'greater', 'standard', 0, 'RD')
  expect_equal(re$n.raw, 3 * n2)
  expect_equal(re$N2.re, ceil_tol(n2))
  expect_equal(re$N1.re, ceiling(2 * ceil_tol(n2)))
})

test_that("coinciding probabilities keep the planned sample size", {
  re <- reestimate(c(0, 0.5), c(0, 0.2), 1, 0.025, 0.8, 'Chisq', 'greater', 'minlike',
                   100L, 0, 'standard', 'group', 30, 30, FALSE, 0, margin.scale = 'RD')
  expect_equal(re$N1.re[1], 30)
  expect_equal(re$N2.re[1], 30)
  expect_equal(re$n.raw[1], 60)
  expect_error(reestimate(0, 0, 1, 0.025, 0.8, 'Chisq', 'greater', 'minlike', 100L, 0,
                          'standard', 'group', NULL, NULL, FALSE, 0,
                          margin.scale = 'RD'), 'coincide')
})

test_that("reestimate keeps the planned size for pairs on the null boundary", {
  # 0.3 - 0.4 + 0.1 = 0 lies on the boundary of the null hypothesis and keeps the planned
  # sizes, whereas two equal probabilities are an ordinary pair under a margin
  re <- reestimate(c(0.3, 0.5), c(0.4, 0.5), 1, 0.025, 0.8, 'Farrington-Manning',
                   'greater', 'minlike', 100L, 0, 'standard', 'nearest', 40, 40, FALSE,
                   0.1, margin.scale = 'RD')
  expect_equal(c(re$N1.re[1], re$N2.re[1]), c(40, 40))
  n2 <- ss_raw_n2(0.5, 0.5, 1, 0.025, 0.8, 'greater', 'standard', 0.1, 'RD')
  expect_equal(re$N2.re[2], floor(n2 + 0.5))
  expect_error(reestimate(0.3, 0.4, 1, 0.025, 0.8, 'Farrington-Manning', 'greater',
                          'minlike', 100L, 0, 'standard', 'nearest', NULL, NULL, FALSE,
                          0.1, margin.scale = 'RD'), 'null boundary')
})

test_that("reestimate passes the search to the exact search", {
  # Reference values from tools/reference/reference_values.py, see test-ss_exact_search.R
  re <- reestimate(c(0.6, 0.6), c(0.3, 0.3), 1, 0.025, 0.85, 'Fisher', 'greater',
                   'minlike', 100L, 0, 'exact', 'group', NULL, NULL, FALSE, 0, 'stable',
                   margin.scale = 'RD')
  expect_equal(c(re$N2.re, re$N2.limit), c(56, 56, 98, 98))
  # The normal approximation ignores the search
  re <- reestimate(0.4, 0.2, 2, 0.025, 0.8, 'Chisq', 'greater', 'minlike', 100L, 0,
                   'standard', 'group', NULL, NULL, FALSE, 0, 'stable',
                   margin.scale = 'RD')
  expect_true(is.na(re$N2.limit))
  expect_equal(re$N2.re, ceil_tol(ss_raw_n2(0.4, 0.2, 2, 0.025, 0.8, 'greater', 'standard',
                                            0, 'RD')))
})

test_that("reestimate keeps the planned size for pairs on a ratio boundary", {
  # 0.24 = 0.8 x 0.3 lies on the boundary p1 = 0.8 p2, and so do two zero rates
  re <- reestimate(c(0.24, 0, 0.5), c(0.3, 0, 0.5), 1, 0.025, 0.8, 'Farrington-Manning',
                   'greater', 'minlike', 100L, 0, 'standard', 'nearest', 40, 40, FALSE,
                   0.8, margin.scale = 'RR')
  expect_equal(re$N2.re[1:2], c(40, 40))
  expect_equal(re$n.raw[1:2], c(80, 80))
  n2 <- ss_raw_n2(0.5, 0.5, 1, 0.025, 0.8, 'greater', 'standard', 0.8, 'RR')
  expect_equal(re$n.raw[3], 2 * n2)
  expect_equal(re$N2.re[3], floor(n2 + 0.5))
  expect_error(reestimate(0, 0, 1, 0.025, 0.8, 'Farrington-Manning', 'greater', 'minlike',
                          100L, 0, 'standard', 'nearest', NULL, NULL, FALSE, 0.8,
                          margin.scale = 'RR'), 'null boundary')
})
