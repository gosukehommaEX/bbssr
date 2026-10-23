test_that("strat_simulate agrees with the exact power of a design with three strata", {
  # Exact reference: enumeration of the interim and second-stage numbers of responders of
  # both groups in every stratum
  set.seed(2026)
  res <- StratPowerBSSR(p = c(0.45, 0.55, 0.35), Delta.A = 0.4, Delta.T = 0.3,
                        N1 = c(4, 4, 4), N2 = c(4, 4, 4), omega = 0.5, r = 1,
                        alpha = 0.1, tar.power = 0.8, Test = 'CMH',
                        ss.method = 'standard', N.max = 36, nsim = 2e5)
  ref <- c(0.547956295988924, 22.297028016)
  expect_lt(abs(res$power.BSSR - ref[1]), 4.5 * res$se.BSSR)
  expect_lt(abs(res$E.N - ref[2]), 0.1)
})

test_that("strat_simulate follows set.seed and returns one value per trial", {
  a <- check_strat_args(c(4, 4, 4), c(4, 4, 4), 0.1, 'CMH', 'greater', 'minlike')
  re.args <- list(Delta.A = 0.4, effect = 'RD', frac = rep(1 / 3, 3), r = c(1, 1, 1),
                  ss.method = 'standard', ss.Test = 'CMH', ss.alpha = 0.1,
                  tar.power = 0.8, N1.plan = c(4L, 4L, 4L), N2.plan = c(4L, 4L, 4L),
                  restricted = FALSE, N.min = NULL, N.max = 36, search = 'crossing',
                  search.limit = c(2, 50))
  set.seed(1)
  s1 <- strat_simulate(c(2L, 2L, 2L), c(2L, 2L, 2L), c(0.6, 0.7, 0.5), c(0.3, 0.4, 0.2),
                       500, 0.1, a, re.args)
  set.seed(1)
  s2 <- strat_simulate(c(2L, 2L, 2L), c(2L, 2L, 2L), c(0.6, 0.7, 0.5), c(0.3, 0.4, 0.2),
                       500, 0.1, a, re.args)
  expect_identical(s1, s2)
  expect_length(s1$reject, 500)
  expect_true(all(s1$N >= 12 & s1$N <= 36))
})
