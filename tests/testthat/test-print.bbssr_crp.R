test_that("print.bbssr_crp reports the outcomes above the level and the largest values", {
  res <- BinaryCondRejectBSSR(Delta.A = 0.3, N1 = 12, N2 = 12, n.interim = c(6, 6), r = 1,
                              alpha = 0.025, tar.power = 0.8, Test = 'Fisher',
                              ss.method = 'standard', theta = c(0.3, 0.5))
  out <- utils::capture.output(ret <- print(res))
  expect_identical(ret, res)
  expect_true(any(grepl('Above the level : 8 for CRP, 0 for CRP.total', out,
                        fixed = TRUE)))
  expect_true(any(grepl('Largest conditional rejection probability', out, fixed = TRUE)))
  expect_true(any(grepl('Type I error rate over 2 values of theta', out, fixed = TRUE)))
  # Without theta the decomposition line is omitted
  res0 <- BinaryCondRejectBSSR(Delta.A = 0.3, N1 = 12, N2 = 12, n.interim = c(6, 6),
                               r = 1, alpha = 0.025, tar.power = 0.8, Test = 'Fisher',
                               ss.method = 'standard')
  out0 <- utils::capture.output(print(res0))
  expect_false(any(grepl('Type I error rate over', out0, fixed = TRUE)))
})
