test_that("PairedBSSR re-estimates the number of pairs of one trial", {
  res <- PairedBSSR(n1 = 10, D1 = 4, Delta.A = 0.25, alpha = 0.05, tar.power = 0.8,
                    Test = 'McNemar', ss.method = 'standard')
  expect_s3_class(res, 'bbssr_pbssr')
  expect_equal(res$N.final, 38)
  expect_equal(res$n.stage2, 28)
  expect_equal(res$psi.hat, 0.4)
  expect_equal(c(res$p10.plan, res$p01.plan), c(0.325, 0.075))
  expect_equal(res$Power, PairedPower(0.325, 0.075, 38, 0.05, 'McNemar')$Power)
  # The final number agrees with the map of PairedPowerBSSR
  map <- attr(PairedPowerBSSR(psi = 0.4, Delta.T = 0.25, Delta.A = 0.25, N = 30,
                              n.interim = 10, alpha = 0.05, tar.power = 0.8,
                              Test = 'McNemar', ss.method = 'standard'), 'reestimation')
  expect_equal(res$N.final, map$N[5])
})

test_that("PairedBSSR reports the conditional power of the rule of Schepers et al.", {
  res <- PairedBSSR(n1 = 10, D1 = 2, Delta.A = 0.25, alpha = 0.05, tar.power = 0.8,
                    Test = 'McNemar-exact', alternative = 'two.sided',
                    ss.method = 'conditional.power', N.max = 45)
  expect_equal(res$N.final, 31)
  expect_gte(res$CP, 0.8)
  rr <- get_paired_rr(30, 0.05, 'McNemar-exact', 'two.sided', 100, 0, FALSE, 0)
  expect_lt(paired_cp(rr, 10L, 2L, 0.25, 1), 0.8)
})

test_that("PairedBSSR stores the full name of the test and warns once", {
  expect_warning(res <- PairedBSSR(n1 = 10, D1 = 4, Delta.A = 0.25, alpha = 0.05,
                                   tar.power = 0.8, Test = 'McNemar-e',
                                   ss.method = 'exact', bb.gamma = 0.001),
                 'only to the McNemar-uncond')
  expect_equal(attr(res, 'Test'), 'McNemar-exact')
  n <- 0
  withCallingHandlers(
    PairedBSSR(n1 = 10, D1 = 4, Delta.A = 0.25, alpha = 0.05, tar.power = 0.8,
               Test = 'McNemar-exact', ss.method = 'exact', bb.gamma = 0.001),
    warning = function(w) {
      n <<- n + 1
      invokeRestart('muffleWarning')
    })
  expect_equal(n, 1)
})

test_that("PairedBSSR applies the restriction and checks the interim data", {
  res <- PairedBSSR(n1 = 10, D1 = 0, Delta.A = 0.25, alpha = 0.05, tar.power = 0.8,
                    Test = 'McNemar', ss.method = 'standard', restricted = TRUE, N = 30)
  expect_equal(res$N.final, 30)
  expect_equal(res$N.re, 23)
  expect_error(PairedBSSR(n1 = 10, D1 = 11, Delta.A = 0.25, alpha = 0.05, tar.power = 0.8,
                          Test = 'McNemar'), 'D1 must be')
  expect_error(PairedBSSR(n1 = 10, D1 = 2, Delta.A = 0.25, alpha = 0.05, tar.power = 0.8,
                          Test = 'McNemar', restricted = TRUE), 'N must be supplied')
  expect_error(PairedBSSR(n1 = 10, D1 = 2, Delta.A = 0.25, alpha = 0.05, tar.power = 0.8,
                          Test = 'McNemar', restricted = TRUE, N = 30.5),
               'single positive integer')
})
