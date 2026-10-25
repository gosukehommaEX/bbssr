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

test_that("print.bbssr_spowerbssr describes a random allocation", {
  set.seed(6)
  res <- StratPowerBSSR(p = c(0.4, 0.6), Delta.A = 0.3, Delta.T = 0.3, N1 = c(6, 6),
                        N2 = c(6, 6), omega = 0.5, r = 1, alpha = 0.05, tar.power = 0.8,
                        Test = 'CMH', ss.method = 'standard', N.max = 30,
                        alloc = 'group', nseq = 5)
  out <- capture.output(print(res))
  expect_true(any(grepl('stratum sizes in the fixed fractions, random group sizes', out,
                        fixed = TRUE)))
  expect_true(any(grepl('first 12 patients, 5 allocation sequences', out, fixed = TRUE)))
  expect_true(any(grepl('se.BSSR', out, fixed = TRUE)))
  expect_false(any(grepl('Interim sizes', out, fixed = TRUE)))
})

test_that("print.bbssr_spowerbssr describes a simulated random allocation", {
  set.seed(23)
  res <- StratPowerBSSR(p = c(0.4, 0.5, 0.6), Delta.A = 0.3, Delta.T = 0.3,
                        N1 = c(2, 2, 1), N2 = c(2, 2, 1), omega = 0.5, r = 1, alpha = 0.05,
                        tar.power = 0.8, Test = 'CMH', ss.method = 'standard',
                        frac = c(0.3, 0.3, 0.4), N.max = 20, alloc = 'stratum',
                        nsim = 500)
  out <- capture.output(print(res))
  expect_true(any(grepl('Interim analysis : first 6 patients', out, fixed = TRUE)))
  expect_false(any(grepl('allocation sequences', out, fixed = TRUE)))
  expect_true(any(grepl('Simulation       : 500 trials per scenario', out, fixed = TRUE)))
  expect_true(any(grepl('se.TRAD', out, fixed = TRUE)))
})
