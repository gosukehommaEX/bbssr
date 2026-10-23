test_that("print.bbssr_ppower reports the power", {
  pw <- PairedPower(c(0.3, 0.35), 0.1, N = 40, alpha = 0.025, Test = 'Paired-Score',
                    margin = 0.05)
  out <- utils::capture.output(ret <- print(pw))
  expect_identical(ret, pw)
  expect_true(any(grepl('Exact power for paired binary data', out)))
  expect_true(any(grepl('Margin      : 0.05', out, fixed = TRUE)))
  expect_true(any(grepl(format(round(pw$Power[2], 4)), out, fixed = TRUE)))
})
