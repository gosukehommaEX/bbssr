test_that("print.bbssr_ppowerbssr reports the design and the power", {
  res <- PairedPowerBSSR(psi = c(0.3, 0.5), Delta.T = 0.25, Delta.A = 0.25, N = 30,
                         n.interim = 10, alpha = 0.05, tar.power = 0.8, Test = 'McNemar',
                         ss.method = 'standard')
  out <- utils::capture.output(ret <- print(res))
  expect_identical(ret, res)
  expect_true(any(grepl('normal approximation (unrestricted)', out, fixed = TRUE)))
  expect_true(any(grepl('N = 30 planned, 10 at the interim analysis', out, fixed = TRUE)))
  expect_true(any(grepl('power.BSSR', out, fixed = TRUE)))
})

test_that("print.bbssr_ppowerbssr prints a subset without its attributes as a data frame", {
  res <- PairedPowerBSSR(psi = c(0.3, 0.5), Delta.T = 0.25, Delta.A = 0.25, N = 30,
                         n.interim = 10, alpha = 0.05, tar.power = 0.8, Test = 'McNemar',
                         ss.method = 'standard')
  # Selecting columns removes the attributes; selecting rows keeps them
  sub <- res[, names(res)]
  expect_null(attr(sub, 'ss.method'))
  out <- utils::capture.output(print(sub))
  expect_false(any(grepl('paired binary data', out, fixed = TRUE)))
  expect_true(any(grepl('power.BSSR', out, fixed = TRUE)))
})
