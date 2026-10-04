# Reference values for Delta.A = 0.3, N1 = N2 = 39, interim sizes 20 and 20, alpha =
# 0.025, target power 0.8, the chi-squared test and re-estimation by the normal
# approximation, computed with an independent Python implementation (numpy and scipy)
# that sums the binomial probabilities over the interim and the second-stage outcomes
tie_args <- list(Delta.A = 0.3, N1 = 39, N2 = 39, n.interim = c(20, 20), r = 1,
                 alpha = 0.025, tar.power = 0.8, Test = 'Chisq', ss.method = 'standard')

test_that("BinaryTypeIErrorBSSR agrees with an independent implementation on a grid", {
  theta <- seq(0.05, 0.95, by = 0.05)
  res <- do.call(BinaryTypeIErrorBSSR, c(tie_args, list(theta = theta, refine = FALSE)))
  expect_s3_class(res, 'bbssr_tie')
  expect_named(res, c('theta', 'TIE.BSSR', 'TIE.TRAD'))
  bssr <- c(0.00587954153053613, 0.0200190561764609, 0.0258863223338634,
            0.0267549931177166, 0.0268501652566501, 0.0265394059647221,
            0.0264742715234961, 0.0270632673880154, 0.0272673261675163,
            0.0271092013571726, 0.0272673261675163, 0.0270632673880154,
            0.0264742715234961, 0.0265394059647222, 0.0268501652566501,
            0.0267549931177166, 0.0258863223338633, 0.0200190561764608,
            0.0058795415305361)
  trad <- c(0.021054107257618, 0.0292884788109491, 0.0266506424787675,
            0.0254801349292422, 0.0260948735515528, 0.0265265508282962,
            0.0255417227643221, 0.0254154435256733, 0.0263080294428053,
            0.0267904681779776, 0.0263080294428052, 0.0254154435256733,
            0.0255417227643221, 0.0265265508282962, 0.0260948735515528,
            0.0254801349292422, 0.0266506424787676, 0.0292884788109491,
            0.021054107257618)
  expect_equal(res$TIE.BSSR, bssr, tolerance = 1e-10)
  expect_equal(res$TIE.TRAD, trad, tolerance = 1e-10)
  m <- attr(res, 'max')
  expect_equal(m$TIE, c(max(bssr), max(trad)), tolerance = 1e-10)
})

test_that("the refined maxima agree with an independent implementation", {
  res <- do.call(BinaryTypeIErrorBSSR, tie_args)
  m <- attr(res, 'max')
  expect_equal(m$TIE, c(0.0272864324696431, 0.0293761105282989), tolerance = 1e-8)
  # Both curves are symmetric about 0.5, so the maximum may be reported at either twin
  expect_equal(pmin(m$theta, 1 - m$theta), c(0.4376846, 0.0935671), tolerance = 1e-4)
  # The refined maximum is at least the largest grid value
  expect_gte(m$TIE[1], max(res$TIE.BSSR))
  expect_gte(m$TIE[2], max(res$TIE.TRAD))
})

test_that("the type I error rate equals the rejection probability of BinaryPowerBSSR", {
  args <- tie_args
  args$n.interim <- NULL
  res <- do.call(BinaryTypeIErrorBSSR, c(args, list(omega = 0.5, theta = c(0.2, 0.5))))
  pw <- BinaryPowerBSSR(p = c(0.2, 0.5), Delta.A = 0.3, Delta.T = 0, N1 = 39, N2 = 39,
                        omega = 0.5, r = 1, alpha = 0.025, tar.power = 0.8, Test = 'Chisq',
                        ss.method = 'standard')
  expect_equal(res$TIE.BSSR, pw$power.BSSR, tolerance = 1e-14)
  expect_equal(res$TIE.TRAD, pw$power.TRAD, tolerance = 1e-14)
})

test_that("BinaryTypeIErrorBSSR validates theta", {
  expect_error(do.call(BinaryTypeIErrorBSSR, c(tie_args, list(theta = c(0.2, 1.2)))),
               'theta')
})
