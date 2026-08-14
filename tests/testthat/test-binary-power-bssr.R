test_that("BinaryPowerBSSR returns a bbssr_powerbssr data frame", {
  res <- BinaryPowerBSSR(
    p = 0.45,
    Delta.A = 0.3, Delta.T = 0.3,
    N1 = 6, N2 = 6, omega = 0.5, r = 1,
    alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
  )
  expect_s3_class(res, 'bbssr_powerbssr')
  expect_named(res, c('p1', 'p2', 'p', 'power.BSSR', 'power.TRAD', 'E.N'))
  expect_equal(nrow(res), 1L)
  expect_true(all(res$power.BSSR >= 0 & res$power.BSSR <= 1))
  expect_true(all(res$power.TRAD >= 0 & res$power.TRAD <= 1))
})

test_that("the arguments of the weighted approach have been removed", {
  base <- list(p = 0.45, Delta.A = 0.3, Delta.T = 0.3, N1 = 6, N2 = 6,
               omega = 0.5, r = 1, alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
  for (dropped in list(list(weighted = TRUE), list(asmd.p1 = 0.6), list(asmd.p2 = 0.3))) {
    expect_error(do.call(BinaryPowerBSSR, c(base, dropped)), 'unused argument',
                 info = names(dropped))
  }
})

test_that("the remaining arguments are all used", {
  # A design is fully determined by Delta.A together with the initial sample sizes,
  # so no argument of the formals may be ignored by the body
  body.text <- paste(deparse(body(BinaryPowerBSSR)), collapse = ' ')
  for (arg in setdiff(names(formals(BinaryPowerBSSR)), '...')) {
    expect_match(body.text, paste0('(?<![\\w.])', arg, '(?![\\w.])'),
                 perl = TRUE, info = arg)
  }
})

test_that("BinaryPowerBSSR is vectorised over the pooled response probability", {
  res <- BinaryPowerBSSR(
    p = seq(0.3, 0.4, by = 0.05),
    Delta.A = 0.3, Delta.T = 0.3,
    N1 = 8, N2 = 8, omega = 0.5, r = 1,
    alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
  )
  expect_equal(nrow(res), 3L)
  expect_equal(res$p1 - res$p2, rep(0.3, 3))
})

test_that("the expected sample size is at least the interim sample size", {
  res <- BinaryPowerBSSR(
    p = 0.35,
    Delta.A = 0.3, Delta.T = 0.3,
    N1 = 10, N2 = 10, omega = 0.5, r = 1,
    alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
  )
  expect_gte(res$E.N, 10)
})

test_that("the restricted rule never enrols fewer patients than the unrestricted rule", {
  args <- list(p = c(0.3, 0.35),
               Delta.A = 0.3, Delta.T = 0.3, N1 = 12, N2 = 12, omega = 0.5, r = 1,
               alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
  unres <- do.call(BinaryPowerBSSR, c(args, list(restricted = FALSE)))
  res <- do.call(BinaryPowerBSSR, c(args, list(restricted = TRUE)))
  expect_true(all(res$E.N >= unres$E.N - 1e-9))
})

test_that("a true treatment effect of zero gives the type I error rate", {
  res <- BinaryPowerBSSR(
    p = c(0.3, 0.5),
    Delta.A = 0.3, Delta.T = 0,
    N1 = 10, N2 = 10, omega = 0.5, r = 1,
    alpha = 0.025, tar.power = 0.8, Test = 'Fisher'
  )
  expect_equal(res$p1, res$p2)
  expect_true(all(res$power.BSSR <= 0.05))
  expect_true(all(res$power.TRAD <= 0.025 + 1e-10))
})

test_that("scenarios outside the unit interval are dropped", {
  res <- BinaryPowerBSSR(
    p = c(0.2, 0.5, 0.95),
    Delta.A = 0.4, Delta.T = 0.4,
    N1 = 6, N2 = 6, omega = 0.5, r = 1,
    alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
  )
  expect_lt(nrow(res), 3L)
  expect_true(all(res$p1 <= 1 & res$p2 >= 0))
})

test_that("BinaryPowerBSSR accepts an interim fraction of one", {
  res <- BinaryPowerBSSR(
    p = 0.45,
    Delta.A = 0.3, Delta.T = 0.3,
    N1 = 6, N2 = 6, omega = 1, r = 1,
    alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
  )
  expect_true(is.finite(res$power.BSSR))
  expect_error(
    BinaryPowerBSSR(p = 0.45,
                    Delta.A = 0.3, Delta.T = 0.3, N1 = 6, N2 = 6, omega = 1.5, r = 1,
                    alpha = 0.025, tar.power = 0.8, Test = 'Chisq'),
    'omega'
  )
})
