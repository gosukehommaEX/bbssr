test_that("print.bbssr_alphaadj reports the adjusted levels", {
  res <- BinaryAlphaAdjBSSR(Delta.A = 0.3, N1 = 12, N2 = 12, omega = 0.5, r = 1,
                            alpha = 0.025, tar.power = 0.8, Test = 'Chisq',
                            ss.method = 'standard', theta = seq(0.1, 0.9, by = 0.1),
                            tol = 1e-4)
  out <- utils::capture.output(ret <- print(res))
  expect_identical(ret, res)
  expect_true(any(grepl('alpha.adj', out)))
  expect_true(any(grepl('final analysis only', out)))
  expect_true(any(grepl('certified over theta in [0.1, 0.9]', out, fixed = TRUE)))
})

test_that("print.bbssr_alphaadj rounds the adjusted levels down", {
  res <- BinaryAlphaAdjBSSR(Delta.A = 0.3, N1 = 12, N2 = 12, omega = 0.5, r = 1,
                            alpha = 0.025, tar.power = 0.8, Test = 'Chisq',
                            ss.method = 'standard', theta = seq(0.1, 0.9, by = 0.1),
                            tol = 1e-4)
  # A level with more digits than printed, which rounding to the nearest would raise
  res$alpha.adj <- c(0.02043759, 0.02077)
  out <- utils::capture.output(print(res))
  expect_true(any(grepl('0.0204375', out, fixed = TRUE)))
  expect_false(any(grepl('0.0204376', out, fixed = TRUE)))
  expect_true(any(grepl('0.02077', out, fixed = TRUE)))
})
