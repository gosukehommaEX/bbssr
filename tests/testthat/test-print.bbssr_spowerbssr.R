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

test_that("print.bbssr_spowerbssr prints a subset without its attributes as a data frame", {
  res <- StratPowerBSSR(p = c(0.4, 0.6), Delta.A = 0.3, Delta.T = 0.3, N1 = c(10, 10),
                        N2 = c(10, 10), omega = 0.5, r = 1, alpha = 0.05, tar.power = 0.8,
                        Test = 'CMH', ss.method = 'standard')
  # Selecting columns removes the attributes; selecting rows keeps them
  for (sub in list(res[, names(res)], res[, c('power.BSSR', 'power.TRAD', 'E.N')])) {
    expect_null(attr(sub, 'ss.method'))
    out <- capture.output(print(sub))
    expect_false(any(grepl('stratified design', out)))
    expect_true(any(grepl('power.BSSR', out, fixed = TRUE)))
  }
})

test_that("print.bbssr_spowerbssr marks the re-estimation from the pooled estimate", {
  res <- StratPowerBSSR(p = c(0.4, 0.6), Delta.A = 0.3, Delta.T = 0.3, N1 = c(10, 10),
                        N2 = c(10, 10), omega = 0.5, r = 1, alpha = 0.05, tar.power = 0.8,
                        Test = 'CMH', ss.method = 'standard', ss.strata = 'pooled')
  out <- capture.output(print(res))
  expect_true(any(grepl('normal approximation (unrestricted, strata pooled)', out,
                        fixed = TRUE)))
})
