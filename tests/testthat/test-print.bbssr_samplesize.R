test_that("print.bbssr_samplesize reports the exact search", {
  out <- utils::capture.output(print(BinarySampleSize(0.6, 0.2, 2, 0.025, 0.9, 'Chisq')))
  expect_true(any(grepl('Exact search     : crossing', out, fixed = TRUE)))
  out <- utils::capture.output(print(BinarySampleSize(0.6, 0.2, 2, 0.025, 0.9, 'Chisq',
                                                      search = 'stable')))
  expect_true(any(grepl('stable up to N2 = 72', out, fixed = TRUE)))
  expect_true(any(grepl('N1 = 46, N2 = 23', out, fixed = TRUE)))
  out <- utils::capture.output(print(BinarySampleSize(0.6, 0.2, 2, 0.025, 0.9, 'Chisq',
                                                      method = 'standard')))
  expect_false(any(grepl('Exact search', out, fixed = TRUE)))
})
