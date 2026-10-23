test_that("print.bbssr_spowerbssr prints the settings and the power", {
  res <- StratPowerBSSR(p = c(0.4, 0.6), Delta.A = 0.3, Delta.T = 0.3, N1 = c(10, 10),
                        N2 = c(10, 10), omega = 0.5, r = 1, alpha = 0.05, tar.power = 0.8,
                        Test = 'CMH', ss.method = 'standard')
  out <- capture.output(r <- print(res))
  expect_identical(r, res)
  expect_true(any(grepl('stratified design', out)))
  expect_true(any(grepl('Planned sizes    : (10, 10) and (10, 10)', out, fixed = TRUE)))
  expect_true(any(grepl('normal approximation (unrestricted)', out, fixed = TRUE)))
})
