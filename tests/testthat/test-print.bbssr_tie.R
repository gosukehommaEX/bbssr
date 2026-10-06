test_that("print.bbssr_tie reports the largest type I error rates", {
  res <- BinaryTypeIErrorBSSR(Delta.A = 0.3, N1 = 12, N2 = 12, omega = 0.5, r = 1,
                              alpha = 0.025, tar.power = 0.8, Test = 'Chisq',
                              ss.method = 'standard', theta = seq(0.1, 0.9, by = 0.1))
  out <- utils::capture.output(ret <- print(res))
  expect_identical(ret, res)
  expect_true(any(grepl('Largest type I error rate', out)))
  expect_true(any(grepl('Fixed sample', out)))
  expect_true(any(grepl('certified over theta in [0.1, 0.9]', out, fixed = TRUE)))
})
