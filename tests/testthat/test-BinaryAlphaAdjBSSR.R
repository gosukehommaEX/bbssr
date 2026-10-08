# Reference values for the design of test-BinaryTypeIErrorBSSR.R, computed with an
# independent Python implementation (numpy and scipy) that bisects over the level of the
# final analysis and maximizes the type I error rate by a grid and a bounded optimization.
# The adjusted levels are reported as the largest values with six significant digits
# below the smallest p-value not rejected, as in report_level()
adj_args <- list(Delta.A = 0.3, N1 = 39, N2 = 39, n.interim = c(20, 20), r = 1,
                 alpha = 0.025, tar.power = 0.8, Test = 'Chisq', ss.method = 'standard')

test_that("BinaryAlphaAdjBSSR agrees with an independent implementation", {
  res <- do.call(BinaryAlphaAdjBSSR, adj_args)
  expect_s3_class(res, 'bbssr_alphaadj')
  expect_named(res, c('Design', 'alpha', 'max.TIE', 'alpha.adj', 'max.TIE.adj', 'theta.adj',
                      'max.TIE.bound', 'max.TIE.adj.bound'))
  expect_equal(res$Design, c('BSSR', 'TRAD'))
  expect_equal(res$max.TIE, c(0.0272864324696431, 0.0293761105282989), tolerance = 1e-8)
  expect_equal(res$alpha.adj, c(0.0204375, 0.02077), tolerance = 1e-12)
  expect_equal(res$max.TIE.adj, c(0.0237418455148, 0.0247660131242), tolerance = 1e-7)
  # The maxima are certified over [0, 1]
  expect_equal(attr(res, 'interval'), c(0, 1))
  expect_true(all(res$max.TIE.bound - res$max.TIE <= 1e-12 + 1e-15))
  expect_true(all(res$max.TIE.adj.bound - res$max.TIE.adj <= 1e-12 + 1e-15))
})

test_that("the adjusted level is certified between the grid points", {
  # With the grid points 0.1 and 0.9 alone the grid misses the largest rates, which the
  # certification over [0.1, 0.9] finds. Reference values from
  # tools/reference/reference_values.py, where every decision is certified
  res <- do.call(BinaryAlphaAdjBSSR, c(adj_args, list(theta = c(0.1, 0.9))))
  expect_equal(attr(res, 'interval'), c(0.1, 0.9))
  expect_equal(res$alpha.adj, c(0.0204375, 0.02077), tolerance = 1e-12)
  expect_equal(res$max.TIE.adj, c(0.0237418455147891, 0.0247660131242263), tolerance = 1e-7)
  # On the grid alone the BSSR design appears to control the type I error rate
  grid <- do.call(BinaryAlphaAdjBSSR,
                  c(adj_args, list(theta = c(0.1, 0.9), maximize = 'grid')))
  expect_equal(grid$alpha.adj[1], 0.025)
  expect_true(all(is.na(grid$max.TIE.bound)))
})

test_that("a failing level is recognised at one grid point without the full evaluation", {
  full <- 0
  real <- refine_max
  local_mocked_bindings(refine_max = function(f, x, y, ...) {
    full <<- full + 1
    real(f, x, y, ...)
  })
  res <- do.call(BinaryAlphaAdjBSSR, adj_args)
  # Each bisection halves an interval of width alpha until it is at most 1e-8 alpha, which
  # takes 27 steps, so without the probe the grid would be evaluated 2 x (1 + 27) times
  expect_gt(full, 2)
  expect_lt(full, 2 * (1 + 27))
  expect_equal(res$alpha.adj, c(0.0204375, 0.02077), tolerance = 1e-12)
})

test_that("the adjusted level controls the type I error rate", {
  res <- do.call(BinaryAlphaAdjBSSR, adj_args)
  tie <- do.call(BinaryTypeIErrorBSSR,
                 c(adj_args[setdiff(names(adj_args), 'alpha')],
                   list(alpha = res$alpha.adj[1], ss.alpha = 0.025)))
  expect_lte(attr(tie, 'max')$TIE[1], 0.025)
  # A level larger by 1e-6 rejects the smallest p-value not rejected and no longer does
  tie.up <- do.call(BinaryTypeIErrorBSSR,
                    c(adj_args[setdiff(names(adj_args), 'alpha')],
                      list(alpha = res$alpha.adj[1] + 1e-6, ss.alpha = 0.025)))
  expect_gt(attr(tie.up, 'max')$TIE[1], 0.025)
})

test_that("a level that already controls the type I error rate is left unchanged", {
  res <- do.call(BinaryAlphaAdjBSSR, modifyList(adj_args, list(Test = 'Fisher')))
  expect_equal(res$alpha.adj, c(0.025, 0.025))
  expect_equal(res$max.TIE.adj, res$max.TIE)
})

test_that("the adjustment of both parts controls the type I error rate", {
  res <- do.call(BinaryAlphaAdjBSSR, c(adj_args, list(adjust = 'both', step = 5e-4)))
  # Reference values from tools/reference/reference_values.py, which certifies the largest
  # rate at 0.025 - k 0.0005 for k = 0, 1, ... and stops at the first level that passes.
  # The levels 0.0245 to 0.0215 fail, and the rate is not monotone in the level
  expect_equal(res$alpha.adj, c(0.021, 0.02077), tolerance = 1e-12)
  expect_equal(res$max.TIE.adj[1], 0.024604289001873, tolerance = 1e-10)
  expect_equal(attr(res, 'adjust'), 'both')
})

test_that("BinaryAlphaAdjBSSR validates tol and step", {
  expect_error(do.call(BinaryAlphaAdjBSSR, c(adj_args, list(tol = 1))), 'tol must be')
  expect_error(do.call(BinaryAlphaAdjBSSR, c(adj_args, list(tol = 0))), 'tol must be')
  expect_error(do.call(BinaryAlphaAdjBSSR, c(adj_args, list(step = 0.025))), 'step must be')
})

test_that("BinaryAlphaAdjBSSR warns when it reports the level 0", {
  # A bisection in which no level passes
  local_mocked_bindings(bisect_level = function(make, assess, certify, alpha, tol) {
    list(level = 0, m0 = list(x = 0.5, y = 0.03, bound = 0.03),
         m = list(x = NA_real_, y = 0, bound = 0))
  })
  expect_warning(res <- BinaryAlphaAdjBSSR(
    Delta.A = 0.3, N1 = 12, N2 = 12, omega = 0.5, r = 1, alpha = 0.025, tar.power = 0.8,
    Test = 'Chisq', ss.method = 'standard', theta = c(0.3, 0.7)
  ), 'BSSR and TRAD design, and the adjusted level is reported as 0')
  expect_equal(res$alpha.adj, c(0, 0))
  expect_true(all(is.na(res$theta.adj)))
})

test_that("BinaryAlphaAdjBSSR passes ref.pvalue to every p-value computation", {
  run <- function(...) {
    BinaryAlphaAdjBSSR(Delta.A = 0.3, N1 = 12, N2 = 12, omega = 0.5, r = 1,
                       alpha = 0.025, tar.power = 0.8, Test = 'Z-pool',
                       ss.method = 'standard', theta = seq(0.05, 0.95, by = 0.05), ...)
  }
  seen <- ref_pvalue_calls(res <- run(ref.pvalue = TRUE))
  expect_gt(nrow(seen), 1)
  expect_true(all(seen$ref.pvalue))
  expect_true(attr(res, 'ref.pvalue'))
  seen <- ref_pvalue_calls(run())
  expect_gt(nrow(seen), 0)
  expect_false(any(seen$ref.pvalue))
})

test_that("BinaryAlphaAdjBSSR evaluates the boundary of a non-inferiority hypothesis", {
  seen <- ref_pvalue_calls(res <- BinaryAlphaAdjBSSR(
    Delta.A = 0, N1 = 54, N2 = 54, n.interim = c(30, 30), r = 1, alpha = 0.025,
    tar.power = 0.8, Test = 'Farrington-Manning', ss.method = 'standard',
    rounding = 'nearest', margin = 0.2, theta = c(0.3, 0.5, 0.7)
  ))
  expect_gt(nrow(seen), 0)
  expect_true(all(seen$margin == 0.2))
  # The largest rates are at least those at the grid points, see
  # test-BinaryTypeIErrorBSSR.R, and the adjusted levels control them
  expect_true(all(res$max.TIE >= c(0.0263276781009502, 0.0258676571158215) - 1e-12))
  expect_true(all(res$alpha.adj < 0.025))
  expect_equal(attr(res, 'margin'), 0.2)
})

test_that("BinaryAlphaAdjBSSR passes tsmethod to every p-value computation", {
  seen <- ref_pvalue_calls(BinaryAlphaAdjBSSR(
    Delta.A = 0.3, N1 = 16, N2 = 8, omega = 0.5, r = 2, alpha = 0.05, tar.power = 0.8,
    Test = 'Fisher', alternative = 'two.sided', tsmethod = 'blaker',
    ss.method = 'standard', theta = seq(0.1, 0.9, by = 0.1)
  ))
  expect_gt(nrow(seen), 1)
  expect_true(all(seen$tsmethod == 'blaker'))
})

test_that("BinaryAlphaAdjBSSR passes search and search.limit to the re-estimation", {
  seen <- search_calls(res <- BinaryAlphaAdjBSSR(
    Delta.A = 0.3, N1 = 8, N2 = 8, n.interim = c(4, 4), r = 1, alpha = 0.025,
    tar.power = 0.8, Test = 'Chisq', search = 'stable', search.limit = c(3, 10),
    theta = c(0.3, 0.5), maximize = 'grid', tol = 1e-4
  ))
  expect_gt(nrow(seen), 0)
  expect_true(all(seen$search == 'stable' & seen$a == 3 & seen$b == 10))
  expect_equal(attr(res, 'search'), 'stable')
})

test_that("BinaryAlphaAdjBSSR evaluates the boundary of a ratio margin", {
  seen <- ref_pvalue_calls(res <- BinaryAlphaAdjBSSR(
    Delta.A = 1, N1 = 100, N2 = 100, n.interim = c(20, 20), r = 1, alpha = 0.025,
    tar.power = 0.8, Test = 'Farrington-Manning', effect = 'RR', ss.method = 'standard',
    rounding = 'nearest', N.max = 300, margin = 0.8, margin.scale = 'RR',
    theta = c(0.3, 0.5, 0.7), maximize = 'grid'
  ))
  expect_gt(nrow(seen), 0)
  expect_true(all(seen$margin.scale == 'RR'))
  # The largest rates on the grid, see test-BinaryTypeIErrorBSSR.R, and the adjusted
  # levels control them
  expect_equal(res$max.TIE, c(0.0260046973465639, 0.025862696947157), tolerance = 1e-10)
  expect_true(all(res$alpha.adj < 0.025))
  expect_true(all(res$max.TIE.adj <= 0.025))
  expect_identical(attr(res, 'margin.scale'), 'RR')
})
