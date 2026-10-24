test_that("print.bbssr_ssamplesize prints the sizes of the strata", {
  res <- StratSampleSize(c(0.5, 0.7), c(0.3, 0.5), c(0.4, 0.6), alpha = 0.025,
                         tar.power = 0.8, Test = 'CMH', method = 'standard')
  out <- capture.output(r <- print(res))
  expect_identical(r, res)
  expect_true(any(grepl('Sample size of a stratified test', out)))
  expect_true(any(grepl('normal approximation', out)))
  expect_true(any(grepl(sprintf('N = %d', res$N), out, fixed = TRUE)))
})

test_that("print.bbssr_ssamplesize describes a random allocation without group sizes", {
  res <- StratSampleSize(c(0.6, 0.8), c(0.2, 0.4), c(0.5, 0.5), alpha = 0.05,
                         tar.power = 0.8, Test = 'CMH', method = 'standard', alloc = 'both')
  out <- capture.output(print(res))
  expect_true(any(grepl('random stratum and group sizes', out, fixed = TRUE)))
  expect_true(any(grepl(sprintf('N = %d', res$N), out, fixed = TRUE)))
  expect_false(any(grepl('N1', out, fixed = TRUE)))
})
