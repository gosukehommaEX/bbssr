test_that("print.bbssr_ssamplesize prints the sizes of the strata", {
  res <- StratSampleSize(c(0.5, 0.7), c(0.3, 0.5), c(0.4, 0.6), alpha = 0.025,
                         tar.power = 0.8, Test = 'CMH', method = 'standard')
  out <- capture.output(r <- print(res))
  expect_identical(r, res)
  expect_true(any(grepl('Sample size of a stratified test', out)))
  expect_true(any(grepl('normal approximation', out)))
  expect_true(any(grepl(sprintf('N = %d', res$N), out, fixed = TRUE)))
})
