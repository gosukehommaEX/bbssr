# Reference values for the odds ratio were computed with an independent Python
# implementation that solves the pooled-rate equation numerically (scipy brentq)
test_that("split_pooled reproduces the risk difference rule of version 2.0.0", {
  p <- c(0.1, 0.35, 0.8)
  sp <- split_pooled(p, 0.3, 2, 'RD')
  expect_identical(sp$p1, p + (1 / 3) * 0.3)
  expect_identical(sp$p2, p - (2 / 3) * 0.3)
})

test_that("split_pooled recovers the risk ratio", {
  sp <- split_pooled(0.3, 1.5, 1, 'RR')
  expect_equal(c(sp$p1, sp$p2), c(0.36, 0.24), tolerance = 1e-12)
  sp <- split_pooled(0.1, 0.4, 1.5, 'RR')
  expect_equal(c(sp$p1, sp$p2), c(0.0625, 0.15625), tolerance = 1e-12)
})

test_that("split_pooled recovers the odds ratio", {
  sp <- split_pooled(0.3, 2, 1, 'OR')
  expect_equal(c(sp$p1, sp$p2), c(0.37171431429143, 0.22828568570857), tolerance = 1e-10)
  sp <- split_pooled(0.42, 0.5, 2, 'OR')
  expect_equal(c(sp$p1, sp$p2), c(0.363439316317354, 0.533121367365293), tolerance = 1e-10)
  sp <- split_pooled(0.8, 3, 2, 'OR')
  expect_equal(c(sp$p1, sp$p2), c(0.862117466393437, 0.675765067213126), tolerance = 1e-10)
})

test_that("the recovered probabilities have the requested pooled value and effect", {
  p <- seq(0.05, 0.95, by = 0.15)
  for (r in c(0.5, 1, 3)) {
    for (psi in c(0.4, 2.5)) {
      sp <- split_pooled(p, psi, r, 'OR')
      expect_equal((r * sp$p1 + sp$p2) / (1 + r), p, tolerance = 1e-12)
      expect_equal(sp$p1 * (1 - sp$p2) / (sp$p2 * (1 - sp$p1)), rep(psi, length(p)),
                   tolerance = 1e-10)
      expect_true(all(sp$p1 > 0 & sp$p1 < 1 & sp$p2 > 0 & sp$p2 < 1))
    }
  }
  # An odds ratio of one leaves the pooled probability unchanged, as do the endpoints
  expect_equal(split_pooled(p, 1, 2, 'OR')$p2, p, tolerance = 1e-14)
  expect_equal(unlist(split_pooled(c(0, 1), 2, 1, 'OR')), c(p11 = 0, p12 = 1, p21 = 0,
                                                            p22 = 1), tolerance = 1e-14)
})
