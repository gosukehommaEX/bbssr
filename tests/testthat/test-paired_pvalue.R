test_that("paired_pvalue obtains the alternative 'less' by exchanging the counts", {
  for (tst in c('McNemar', 'McNemar-exact', 'McNemar-midP', 'McNemar-uncond')) {
    expect_equal(paired_pvalue(9, tst, 'less', 30, 0, FALSE, 0),
                 t(paired_pvalue(9, tst, 'greater', 30, 0, FALSE, 0)), info = tst)
  }
  expect_equal(paired_pvalue(9, 'Paired-Score', 'less', 30, 0, FALSE, 0.1),
               t(paired_pvalue(9, 'Paired-Score', 'greater', 30, 0, FALSE, 0.1)))
})

test_that("paired_pvalue takes the larger one-sided p-value for equivalence", {
  P <- paired_pvalue(12, 'Paired-Score', 'equivalence', 30, 0, FALSE, c(0.1, 0.15))
  lo <- paired_pvalue(12, 'Paired-Score', 'greater', 30, 0, FALSE, 0.1)
  up <- paired_pvalue(12, 'Paired-Score', 'less', 30, 0, FALSE, 0.15)
  expect_equal(P, pmax(lo, up))
  cells <- cbind(c(5, 6, 3, 0, 2, 4), c(1, 2, 0, 0, 3, 4)) + 1
  expect_equal(P[cells], c(0.829587827230053, 0.788380835652461, 0.834012266458632,
                           0.124106539494962, 0.464034985075496, 0.335099118237478),
               tolerance = 1e-10)
})

test_that("paired_pvalue of the asymptotic tests refers Z to the normal distribution", {
  Z <- paired_score_z(10, 0)
  expect_equal(paired_pvalue(10, 'McNemar', 'greater', 30, 0, FALSE, 0),
               stats::pnorm(Z, lower.tail = FALSE))
  expect_equal(paired_pvalue(10, 'McNemar', 'two.sided', 30, 0, FALSE, 0),
               pmin(2 * stats::pnorm(abs(Z), lower.tail = FALSE), 1))
  # The score test with margin 0 is the McNemar test
  expect_equal(paired_pvalue(10, 'Paired-Score', 'greater', 30, 0, FALSE, 0),
               paired_pvalue(10, 'McNemar', 'greater', 30, 0, FALSE, 0))
})
