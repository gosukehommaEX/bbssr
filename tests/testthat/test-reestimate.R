test_that("reestimate evaluates each distinct pair once and reproduces BinarySampleSize", {
  hat.p1 <- c(0.5, 0.6, 0.5)
  hat.p2 <- c(0.2, 0.3, 0.2)
  re <- reestimate(hat.p1, hat.p2, 1, 0.025, 0.8, 'Fisher', 'greater', 'minlike', 100L,
                   0, 'exact', 'group', NULL, NULL, FALSE)
  expect_named(re, c('n.raw', 'N1.re', 'N2.re'))
  expect_true(all(is.na(re$n.raw)))
  for (i in 1:2) {
    ss <- BinarySampleSize(hat.p1[i], hat.p2[i], 1, 0.025, 0.8, 'Fisher')
    expect_equal(c(re$N1.re[i], re$N2.re[i]), c(ss$N1, ss$N2))
  }
  expect_equal(re[3, ], re[1, ], ignore_attr = TRUE)
})

test_that("the normal methods return the unrounded total", {
  re <- reestimate(0.4, 0.2, 2, 0.025, 0.8, 'Chisq', 'greater', 'minlike', 100L, 0,
                   'standard', 'group', NULL, NULL, FALSE)
  n2 <- ss_raw_n2(0.4, 0.2, 2, 0.025, 0.8, 'greater', 'standard')
  expect_equal(re$n.raw, 3 * n2)
  expect_equal(re$N2.re, ceil_tol(n2))
  expect_equal(re$N1.re, ceiling(2 * ceil_tol(n2)))
})

test_that("coinciding probabilities keep the planned sample size", {
  re <- reestimate(c(0, 0.5), c(0, 0.2), 1, 0.025, 0.8, 'Chisq', 'greater', 'minlike',
                   100L, 0, 'standard', 'group', 30, 30, FALSE)
  expect_equal(re$N1.re[1], 30)
  expect_equal(re$N2.re[1], 30)
  expect_equal(re$n.raw[1], 60)
  expect_error(reestimate(0, 0, 1, 0.025, 0.8, 'Chisq', 'greater', 'minlike', 100L, 0,
                          'standard', 'group', NULL, NULL, FALSE), 'coincide')
})
