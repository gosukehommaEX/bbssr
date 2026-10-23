test_that("print.bbssr_stie prints the largest rates and how they were located", {
  tie <- StratTypeIErrorBSSR(Delta.A = 0.3, N1 = c(6, 6), N2 = c(6, 6), omega = 0.5, r = 1,
                             alpha = 0.05, tar.power = 0.8, Test = 'CMH',
                             ss.method = 'standard', theta = c(0.2, 0.5, 0.8))
  out <- capture.output(r <- print(tie))
  expect_identical(r, tie)
  expect_true(any(grepl('certified over (theta.S1, theta.S2) in [0.2, 0.8] x [0.2, 0.8]',
                        out, fixed = TRUE)))
  expect_true(any(grepl('Fixed sample', out)))
  expect_true(any(grepl('bound', out)))
})
