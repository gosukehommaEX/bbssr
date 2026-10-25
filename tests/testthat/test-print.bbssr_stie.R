test_that("print.bbssr_stie prints the largest rates and how they were located", {
  tie <- StratTypeIErrorBSSR(Delta.A = 0.3, N1 = c(6, 6), N2 = c(6, 6), omega = 0.5, r = 1,
                             alpha = 0.05, tar.power = 0.8, Test = 'CMH',
                             ss.method = 'standard', theta = c(0.2, 0.5, 0.8))
  out <- capture.output(r <- print(tie))
  expect_identical(r, tie)
  expect_true(any(grepl('certified over (theta.S1, theta.S2) in [0.2, 0.8] x [0.2, 0.8]',
                        out, fixed = TRUE)))
  expect_true(any(grepl('Fixed sample', out)))
  expect_true(any(grepl('bound', out)))
})

test_that("print.bbssr_stie prints a subset without its attributes as a data frame", {
  tie <- StratTypeIErrorBSSR(Delta.A = 0.3, N1 = c(6, 6), N2 = c(6, 6), omega = 0.5, r = 1,
                             alpha = 0.05, tar.power = 0.8, Test = 'CMH',
                             ss.method = 'standard', theta = c(0.2, 0.5, 0.8))
  # Selecting columns removes the attributes; selecting rows keeps them
  sub <- tie[, names(tie)]
  expect_null(attr(sub, 'max'))
  out <- capture.output(print(sub))
  expect_false(any(grepl('Type I error rate of blinded', out)))
  expect_true(any(grepl('TIE.BSSR', out, fixed = TRUE)))
})

test_that("print.bbssr_stie describes a random allocation and the standard error", {
  set.seed(7)
  res <- StratTypeIErrorBSSR(Delta.A = 0.3, N1 = c(6, 6), N2 = c(6, 6), omega = 0.5, r = 1,
                             alpha = 0.05, tar.power = 0.8, Test = 'CMH',
                             ss.method = 'standard', N.max = 30, theta = c(0.3, 0.7),
                             alloc = 'stratum', nseq = 5)
  out <- capture.output(print(res))
  expect_true(any(grepl('permuted blocks of 4 within the strata', out, fixed = TRUE)))
  expect_true(any(grepl('Sequences        : 5', out, fixed = TRUE)))
  expect_true(any(grepl(' se', out, fixed = TRUE)))
})

test_that("print.bbssr_stie describes a simulated random allocation", {
  set.seed(24)
  res <- StratTypeIErrorBSSR(Delta.A = 0.3, N1 = c(2, 2, 1), N2 = c(2, 2, 1), omega = 0.5,
                             r = 1, alpha = 0.05, tar.power = 0.8, Test = 'CMH',
                             ss.method = 'standard', frac = c(0.3, 0.3, 0.4), N.max = 20,
                             theta = c(0.3, 0.6), alloc = 'both', nsim = 500)
  out <- capture.output(print(res))
  expect_true(any(grepl('random stratum and group sizes', out, fixed = TRUE)))
  expect_true(any(grepl('Interim analysis : first 6 patients', out, fixed = TRUE)))
  expect_false(any(grepl('Sequences', out, fixed = TRUE)))
})
