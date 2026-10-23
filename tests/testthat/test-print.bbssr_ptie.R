test_that("print.bbssr_ptie reports the largest type I error rates", {
  tie <- PairedTypeIErrorBSSR(Delta.A = 0.25, N = 30, n.interim = 10, alpha = 0.05,
                              tar.power = 0.8, Test = 'McNemar', ss.method = 'standard',
                              psi = seq(0.1, 0.9, by = 0.1))
  out <- utils::capture.output(ret <- print(tie))
  expect_identical(ret, tie)
  expect_true(any(grepl('certified over psi in [0.1, 0.9]', out, fixed = TRUE)))
  expect_true(any(grepl('9 values of psi in [0.1, 0.9]', out, fixed = TRUE)))
  expect_true(any(grepl('Fixed sample', out, fixed = TRUE)))
})
