# Reference values for Delta.A = 0.3, N1 = N2 = 39, interim sizes 20 and 20, alpha =
# 0.025, target power 0.8, the chi-squared test and re-estimation by the normal
# approximation, computed with an independent Python implementation (numpy and scipy)
# that sums the binomial probabilities over the interim and the second-stage outcomes
tie_args <- list(Delta.A = 0.3, N1 = 39, N2 = 39, n.interim = c(20, 20), r = 1,
                 alpha = 0.025, tar.power = 0.8, Test = 'Chisq', ss.method = 'standard')

test_that("BinaryTypeIErrorBSSR agrees with an independent implementation on a grid", {
  theta <- seq(0.05, 0.95, by = 0.05)
  res <- do.call(BinaryTypeIErrorBSSR, c(tie_args, list(theta = theta, maximize = 'grid')))
  expect_s3_class(res, 'bbssr_tie')
  expect_named(res, c('theta', 'p1', 'p2', 'TIE.BSSR', 'TIE.TRAD'))
  # Without a margin both groups share the response probability theta
  expect_identical(res$p1, res$theta)
  expect_identical(res$p2, res$theta)
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
  expect_true(all(is.na(m$bound)))
  expect_equal(attr(res, 'maximize'), 'grid')
  expect_null(attr(res, 'interval'))
})

test_that("the certified maxima agree with an independent implementation", {
  res <- do.call(BinaryTypeIErrorBSSR, tie_args)
  m <- attr(res, 'max')
  expect_equal(attr(res, 'maximize'), 'certified')
  expect_equal(attr(res, 'interval'), c(0, 1))
  expect_equal(m$TIE, c(0.0272864324696176, 0.0293761105282982), tolerance = 1e-10)
  expect_true(all(m$bound - m$TIE <= 1e-12 + 1e-15))
  # Both curves are symmetric about 0.5, so the maximum may be reported at either twin
  expect_equal(pmin(m$theta, 1 - m$theta), c(0.4376846, 0.0935671), tolerance = 1e-4)
  # The refinement of the largest local maxima on the grid finds the same values here
  ref <- attr(do.call(BinaryTypeIErrorBSSR, c(tie_args, list(maximize = 'refined'))), 'max')
  expect_equal(ref$TIE, c(0.0272864324696431, 0.0293761105282989), tolerance = 1e-8)
  expect_true(all(is.na(ref$bound)))
  # The refined maxima, computed by summation, do not exceed the bounds
  expect_true(all(ref$TIE <= m$bound + 1e-15))
})

test_that("the certified maximum covers the interval between the grid points", {
  res <- do.call(BinaryTypeIErrorBSSR, c(tie_args, list(theta = c(0.15, 0.85))))
  m <- attr(res, 'max')
  expect_equal(attr(res, 'interval'), c(0.15, 0.85))
  expect_equal(m$TIE, c(0.0272864324696353, 0.0267904681779775), tolerance = 1e-10)
  expect_true(all(m$bound - m$TIE <= 1e-12 + 1e-15))
  # Both maxima lie between the two grid points, where the grid cannot see them
  expect_true(all(m$theta > 0.15 & m$theta < 0.85))
  expect_gt(m$TIE[1], max(res$TIE.BSSR))
  expect_gt(m$TIE[2], max(res$TIE.TRAD))
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

test_that("ref.pvalue reaches every rejection region and removes the grid excess", {
  run <- function(...) {
    BinaryTypeIErrorBSSR(Delta.A = 0.3, N1 = 32, N2 = 32, omega = 0.5, r = 1,
                         alpha = 0.025, tar.power = 0.8, Test = 'Z-pool',
                         ss.method = 'standard', ...)
  }
  seen <- ref_pvalue_calls(ref <- run(ref.pvalue = TRUE))
  expect_gt(nrow(seen), 1)
  expect_true(all(seen$ref.pvalue))
  expect_true(attr(ref, 'ref.pvalue'))
  grid <- run()
  # Largest type I error rate of the fixed design with 32 patients per group, certified
  # by tools/reference/reference_values.py
  expect_equal(attr(grid, 'max')$TIE[2], 0.0250058337172828, tolerance = 1e-9)
  expect_equal(attr(ref, 'max')$TIE[2], 0.023344357650963, tolerance = 1e-9)
})

test_that("BinaryTypeIErrorBSSR evaluates the boundary of a non-inferiority hypothesis", {
  seen <- ref_pvalue_calls(tie <- BinaryTypeIErrorBSSR(
    Delta.A = 0, N1 = 54, N2 = 54, n.interim = c(30, 30), r = 1, alpha = 0.025,
    tar.power = 0.8, Test = 'Farrington-Manning', ss.method = 'standard',
    rounding = 'nearest', margin = 0.2, theta = c(0.05, 0.3, 0.5, 0.7, 0.95),
    maximize = 'grid'
  ))
  expect_gt(nrow(seen), 0)
  expect_true(all(seen$margin == 0.2))
  # theta = 0.05 and 0.95 put a response probability outside the unit interval
  expect_equal(tie$theta, c(0.3, 0.5, 0.7))
  expect_equal(tie$p1 - tie$p2, rep(-0.2, 3))
  # Reference values from tools/reference/reference_values.py, BSSR and fixed design in
  # turn for each value of theta
  expect_equal(c(rbind(tie$TIE.BSSR, tie$TIE.TRAD)),
               c(0.0263276781009502, 0.0258676571158215, 0.0245342677114993,
                 0.0229853663986697, 0.0263276781009502, 0.0258676571158216),
               tolerance = 1e-10)
  expect_equal(attr(tie, 'margin'), 0.2)
})

test_that("the largest type I error rate on a non-inferiority boundary is certified", {
  tie <- BinaryTypeIErrorBSSR(
    Delta.A = 0, N1 = 54, N2 = 54, n.interim = c(30, 30), r = 1, alpha = 0.025,
    tar.power = 0.8, Test = 'Farrington-Manning', ss.method = 'standard',
    rounding = 'nearest', margin = 0.2
  )
  # The default grid spans [0, 1], and the boundary p1 - p2 = -0.2 has theta in [0.1, 0.9]
  expect_equal(attr(tie, 'interval'), c(0.1, 0.9))
  expect_equal(range(tie$theta), c(0.1, 0.9))
  m <- attr(tie, 'max')
  # Reference values from tools/reference/reference_values.py
  expect_equal(m$TIE, c(0.0263983576247002, 0.0281499487826859), tolerance = 1e-10)
  expect_true(all(m$bound - m$TIE <= 1e-12 + 1e-15))
})

test_that("the certified maximum agrees with a dense grid for unequal groups and 'less'", {
  # Groups of different sizes and the boundary p1 - p2 = 0.3 of the hypothesis of
  # alternative = 'less', on which p1 = theta + 0.1 and p2 = theta - 0.2, so theta runs
  # from 0.2 to 0.9. The certification uses the Bernstein coefficients and the refinement
  # on a grid of step 0.001 sums the binomial probabilities, two separate computations
  args <- list(Delta.A = 0, N1 = 24, N2 = 12, n.interim = c(12, 6), r = 2, alpha = 0.025,
               tar.power = 0.8, Test = 'Farrington-Manning', alternative = 'less',
               ss.method = 'standard', margin = 0.3)
  cert <- do.call(BinaryTypeIErrorBSSR, c(args, list(theta = c(0.2, 0.9))))
  expect_equal(attr(cert, 'interval'), c(0.2, 0.9), tolerance = 1e-14)
  expect_equal(cert$p1 - cert$p2, c(0.3, 0.3), tolerance = 1e-14)
  dense <- do.call(BinaryTypeIErrorBSSR,
                   c(args, list(theta = seq(0.2, 0.9, by = 0.001), maximize = 'refined')))
  mc <- attr(cert, 'max')
  md <- attr(dense, 'max')
  expect_gt(min(md$TIE), 0.001)
  expect_equal(mc$TIE, md$TIE, tolerance = 1e-10)
  expect_true(all(md$TIE <= mc$bound + 1e-15))
  expect_true(all(mc$bound - mc$TIE <= 1e-12 + 1e-15))
})

test_that("BinaryTypeIErrorBSSR passes tsmethod to every rejection region", {
  seen <- ref_pvalue_calls(BinaryTypeIErrorBSSR(
    Delta.A = 0.3, N1 = 16, N2 = 8, omega = 0.5, r = 2, alpha = 0.05, tar.power = 0.8,
    Test = 'Fisher', alternative = 'two.sided', tsmethod = 'blaker',
    ss.method = 'standard', theta = c(0.3, 0.5)
  ))
  expect_gt(nrow(seen), 1)
  expect_true(all(seen$tsmethod == 'blaker'))
})

test_that("BinaryTypeIErrorBSSR passes search and search.limit to the re-estimation", {
  seen <- search_calls(res <- BinaryTypeIErrorBSSR(
    Delta.A = 0.3, N1 = 8, N2 = 8, n.interim = c(4, 4), r = 1, alpha = 0.025,
    tar.power = 0.8, Test = 'Chisq', search = 'stable', search.limit = c(3, 10),
    theta = c(0.3, 0.5), maximize = 'grid'
  ))
  expect_gt(nrow(seen), 0)
  expect_true(all(seen$search == 'stable' & seen$a == 3 & seen$b == 10))
  expect_equal(attr(res, 'search'), 'stable')
  expect_false(anyNA(attr(res, 'reestimation')$N2.limit))
})

test_that("BinaryTypeIErrorBSSR evaluates the boundary of a ratio margin", {
  args <- list(Delta.A = 1, N1 = 100, N2 = 100, n.interim = c(20, 20), r = 1,
               alpha = 0.025, tar.power = 0.8, Test = 'Farrington-Manning', effect = 'RR',
               ss.method = 'standard', rounding = 'nearest', N.max = 300, margin = 0.8,
               margin.scale = 'RR')
  seen <- ref_pvalue_calls(tie <- do.call(BinaryTypeIErrorBSSR, c(args, list(
    theta = c(0.3, 0.5, 0.7, 0.95), maximize = 'grid'
  ))))
  expect_true(all(seen$margin.scale == 'RR'))
  # theta = 0.95 puts the response probability of group 2 above 1
  expect_equal(tie$theta, c(0.3, 0.5, 0.7))
  expect_equal(tie$p1, 0.8 * tie$p2)
  # Reference values from tools/reference/reference_values.py, BSSR and fixed design in
  # turn for each value of theta
  expect_equal(c(rbind(tie$TIE.BSSR, tie$TIE.TRAD)),
               c(0.0255030250799468, 0.025862696947157, 0.0253098748918834,
                 0.025463971862289, 0.0260046973465639, 0.0249123070548097),
               tolerance = 1e-10)
  expect_identical(attr(tie, 'margin.scale'), 'RR')
  # Certified over the whole boundary, on which theta runs from 0 to 0.9 and p2 from 0 to
  # 1. The largest rate of the re-estimation design lies at the end, where few patients
  # are recruited
  cert <- do.call(BinaryTypeIErrorBSSR, args)
  expect_equal(attr(cert, 'interval'), c(0, 0.9), tolerance = 1e-14)
  m <- attr(cert, 'max')
  expect_equal(m$TIE, c(0.0370092046191537, 0.0262750446720697), tolerance = 1e-10)
  expect_equal(m$theta[1], 0.9, tolerance = 1e-9)
  expect_true(all(m$bound - m$TIE <= 1e-12 + 1e-15))
})

test_that("the certified maximum on a ratio boundary agrees with a dense grid", {
  # Groups of different sizes and the boundary p1 = 1.5 p2 of the hypothesis of
  # alternative = 'less', on which p2 = 3 theta / 4 and p1 = 1.125 theta, so theta can run
  # up to 8 / 9. The certification over [0.05, 0.85] uses the Bernstein coefficients and
  # the refinement on a grid of step 0.001 sums the binomial probabilities, two separate
  # computations
  args <- list(Delta.A = 1, N1 = 24, N2 = 12, n.interim = c(12, 6), r = 2, alpha = 0.025,
               tar.power = 0.8, Test = 'Farrington-Manning', alternative = 'less',
               effect = 'RR', ss.method = 'standard', N.max = 120, margin = 1.5,
               margin.scale = 'RR')
  cert <- do.call(BinaryTypeIErrorBSSR, c(args, list(theta = c(0.05, 0.85))))
  expect_equal(attr(cert, 'interval'), c(0.05, 0.85), tolerance = 1e-14)
  expect_equal(cert$p1, 1.5 * cert$p2, tolerance = 1e-14)
  dense <- do.call(BinaryTypeIErrorBSSR,
                   c(args, list(theta = seq(0.05, 0.85, by = 0.001), maximize = 'refined')))
  mc <- attr(cert, 'max')
  md <- attr(dense, 'max')
  expect_gt(min(md$TIE), 0.001)
  expect_equal(mc$TIE, md$TIE, tolerance = 1e-10)
  expect_true(all(md$TIE <= mc$bound + 1e-15))
  expect_true(all(mc$bound - mc$TIE <= 1e-12 + 1e-15))
})
