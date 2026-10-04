# Reference values for p = 0.4 and Delta.T = 0.3, computed with an independent Python
# implementation (numpy and scipy) from the convolution of the two interim binomial
# distributions
test_that("summary reports the distribution of the final sample size", {
  res <- BinaryPowerBSSR(p = 0.4, Delta.A = 0.3, Delta.T = 0.3, N1 = 39, N2 = 39,
                         n.interim = c(20, 20), r = 1, alpha = 0.025, tar.power = 0.8,
                         Test = 'Chisq', ss.method = 'standard')
  s <- summary(res)
  expect_named(s, c('p1', 'p2', 'p', 'E.N', 'SD.N', 'N.q25', 'N.q50', 'N.q75', 'P.N.max'))
  expect_equal(s$E.N, 80.3764079043752, tolerance = 1e-10)
  expect_equal(s$E.N, res$E.N, tolerance = 1e-12)
  expect_equal(s$SD.N, 5.59614424051928, tolerance = 1e-10)
  expect_equal(c(s$N.q25, s$N.q50, s$N.q75), c(78, 82, 84))
  expect_gte(s$P.N.max, 0)
})

test_that("the distribution of the final sample size sums to one for every scenario", {
  res <- BinaryPowerBSSR(p = c(0.3, 0.3, 0.5), Delta.A = 0.3, Delta.T = 0.3, N1 = 12,
                         N2 = 12, omega = 0.5, r = 1, alpha = 0.025, tar.power = 0.8,
                         Test = 'Chisq')
  d <- attr(res, 'N.dist')
  expect_equal(as.vector(tapply(d$prob, d$scenario, sum)), rep(1, 3), tolerance = 1e-12)
  s <- summary(res, probs = c(0.1, 0.9))
  expect_named(s, c('p1', 'p2', 'p', 'E.N', 'SD.N', 'N.q10', 'N.q90', 'P.N.max'))
  expect_equal(s[1, ], s[2, ], ignore_attr = TRUE)
  expect_error(summary(res, probs = 0), 'probs')
})
