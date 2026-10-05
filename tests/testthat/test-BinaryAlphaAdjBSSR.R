# Reference values for the design of test-BinaryTypeIErrorBSSR.R, computed with an
# independent Python implementation (numpy and scipy) that bisects over the level of the
# final analysis and maximizes the type I error rate by a grid and a bounded optimization
adj_args <- list(Delta.A = 0.3, N1 = 39, N2 = 39, n.interim = c(20, 20), r = 1,
                 alpha = 0.025, tar.power = 0.8, Test = 'Chisq', ss.method = 'standard')

test_that("BinaryAlphaAdjBSSR agrees with an independent implementation", {
  res <- do.call(BinaryAlphaAdjBSSR, adj_args)
  expect_s3_class(res, 'bbssr_alphaadj')
  expect_equal(res$Design, c('BSSR', 'TRAD'))
  expect_equal(res$max.TIE, c(0.0272864324696431, 0.0293761105282989), tolerance = 1e-8)
  expect_equal(res$alpha.adj, c(0.0204375745496, 0.0207700483967), tolerance = 1e-7)
  expect_equal(res$max.TIE.adj, c(0.0237418455148, 0.0247660131242), tolerance = 1e-7)
  expect_true(all(res$max.TIE.adj <= 0.025))
})

test_that("the adjusted level controls the type I error rate", {
  res <- do.call(BinaryAlphaAdjBSSR, adj_args)
  tie <- do.call(BinaryTypeIErrorBSSR,
                 c(adj_args[setdiff(names(adj_args), 'alpha')],
                   list(alpha = res$alpha.adj[1], ss.alpha = 0.025)))
  expect_lte(attr(tie, 'max')$TIE[1], 0.025)
  # A slightly larger level no longer does
  tie.up <- do.call(BinaryTypeIErrorBSSR,
                    c(adj_args[setdiff(names(adj_args), 'alpha')],
                      list(alpha = res$alpha.adj[1] * (1 + 1e-6), ss.alpha = 0.025)))
  expect_gt(attr(tie.up, 'max')$TIE[1], 0.025)
})

test_that("a level that already controls the type I error rate is left unchanged", {
  res <- do.call(BinaryAlphaAdjBSSR, modifyList(adj_args, list(Test = 'Fisher')))
  expect_equal(res$alpha.adj, c(0.025, 0.025))
  expect_equal(res$max.TIE.adj, res$max.TIE)
})

test_that("the adjustment of both parts controls the type I error rate", {
  res <- do.call(BinaryAlphaAdjBSSR, c(adj_args, list(adjust = 'both', step = 5e-4)))
  expect_lte(res$max.TIE.adj[1], 0.025)
  expect_lt(res$alpha.adj[1], 0.025)
  expect_equal(attr(res, 'adjust'), 'both')
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
  expect_false(any(seen$ref.pvalue))
})
