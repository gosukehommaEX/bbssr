test_that("paired_pvalue obtains the alternative 'less' by exchanging the counts", {
  for (tst in c('McNemar', 'McNemar-cc', 'McNemar-exact', 'McNemar-midP',
                'McNemar-uncond')) {
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

test_that("paired_pvalue of the corrected McNemar test agrees with stats::mcnemar.test", {
  N <- 15
  P <- paired_pvalue(N, 'McNemar-cc', 'two.sided', 30, 0, FALSE, 0)
  ref <- matrix(NA_real_, N + 1, N + 1)
  for (a in 0:N) {
    for (b in 0:(N - a)) {
      # stats::mcnemar.test() gives NaN without discordant pairs
      if (a + b > 0) {
        ref[a + 1, b + 1] <- stats::mcnemar.test(matrix(c(0, a, b, 0), 2),
                                                 correct = TRUE)$p.value
      }
    }
  }
  ok <- !is.na(ref)
  expect_equal(sum(ok), (N + 1) * (N + 2) / 2 - 1)
  expect_equal(P[ok], ref[ok], tolerance = 1e-10)
  expect_equal(P[1, 1], 1)
  expect_true(all(is.na(P[outer(0:N, 0:N, '+') > N])))
  # The two-sided p-value is twice the smaller one-sided p-value, at most 1, and the
  # one-sided p-value refers (a - b - 1) / sqrt(a + b) to the normal distribution
  G <- paired_pvalue(N, 'McNemar-cc', 'greater', 30, 0, FALSE, 0)
  L <- paired_pvalue(N, 'McNemar-cc', 'less', 30, 0, FALSE, 0)
  expect_equal(P, pmin(2 * pmin(G, L), 1))
  expect_equal(G[9, 3], stats::pnorm((8 - 2 - 1) / sqrt(10), lower.tail = FALSE))
  expect_equal(G[1, 1], 1)
})
