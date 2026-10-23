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

test_that("print shows both boundaries of a test of equivalence", {
  tie <- BinaryTypeIErrorBSSR(Delta.A = 0, N1 = 30, N2 = 30, n.interim = c(15, 15),
                              r = 1, alpha = 0.05, tar.power = 0.8,
                              Test = 'Farrington-Manning', alternative = 'equivalence',
                              ss.method = 'standard', margin = 0.3,
                              theta = c(0.3, 0.5, 0.7), maximize = 'grid')
  out <- capture.output(print(tie))
  expect_true(any(grepl('3 values of theta on the lower boundary and 3 on the upper', out)))
  expect_true(any(grepl('Boundary', out)))
  expect_true(any(grepl('-0.3 to 0.3 (risk difference)', out, fixed = TRUE)))
})
