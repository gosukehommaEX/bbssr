test_that("PairedAlphaAdjBSSR agrees with the reference adjusted levels", {
  adj <- PairedAlphaAdjBSSR(Delta.A = 0.25, N = 30, n.interim = 10, alpha = 0.05,
                            tar.power = 0.8, Test = 'McNemar', ss.method = 'standard')
  expect_s3_class(adj, 'bbssr_alphaadj')
  expect_equal(adj$max.TIE, c(0.0525269607489134, 0.0558117064555376), tolerance = 1e-9)
  expect_equal(adj$alpha.adj, c(0.0421657, 0.0416322), tolerance = 1e-12)
  expect_equal(adj$max.TIE.adj, c(0.049970120597453, 0.0493685733526942), tolerance = 1e-8)
  expect_true(all(adj$max.TIE.adj.bound <= 0.05))
  expect_true(all(adj$psi.adj >= 0 & adj$psi.adj <= 1))
  out <- utils::capture.output(print(adj))
  expect_true(any(grepl('certified over psi in [0, 1]', out, fixed = TRUE)))
})

test_that("PairedAlphaAdjBSSR keeps the level of the exact conditional test", {
  adj <- PairedAlphaAdjBSSR(Delta.A = 0.25, N = 30, n.interim = 10, alpha = 0.05,
                            tar.power = 0.8, Test = 'McNemar-exact',
                            alternative = 'two.sided', ss.method = 'standard',
                            psi = seq(0, 1, by = 0.05))
  expect_equal(adj$alpha.adj, c(0.05, 0.05))
})

test_that("PairedAlphaAdjBSSR checks its arguments", {
  run <- function(...) {
    PairedAlphaAdjBSSR(Delta.A = 0.25, N = 30, n.interim = 10, alpha = 0.05,
                       tar.power = 0.8, Test = 'McNemar', ss.method = 'standard', ...)
  }
  expect_error(run(tol = 1), 'tol must be')
  expect_error(run(step = 0.06), 'step must be')
  expect_error(run(psi = -0.1), 'psi must be')
})
