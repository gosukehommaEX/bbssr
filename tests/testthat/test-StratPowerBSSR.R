pb <- function(...) {
  args <- list(p = rbind(c(0.4, 0.6), c(0.4, 0.6)), Delta.A = 0.3, Delta.T = c(0.3, 0),
               N1 = c(10, 10), N2 = c(10, 10), omega = 0.5, r = 1, alpha = 0.05,
               tar.power = 0.8, Test = 'CMH', ss.method = 'standard')
  new <- list(...)
  args[names(new)] <- new
  do.call(StratPowerBSSR, args)
}

test_that("StratPowerBSSR agrees with the reference values with two strata", {
  # Scenario 1: p1 = (0.55, 0.75), p2 = (0.25, 0.45); scenario 2: no effect
  one <- pb(p = c(0.4, 0.6), Delta.T = 0.3)
  two <- pb(p = c(0.4, 0.6), Delta.T = 0)
  expect_equal(c(one$power.BSSR, one$E.N, two$power.BSSR, two$E.N),
               c(0.770441457550367, 59.7402124907309, 0.0502722588565529, 59.1902574387837),
               tolerance = 1e-10)
  ex <- function(dT) {
    pb(p = c(0.4, 0.6), Delta.A = 0.4, Delta.T = dT, N1 = c(8, 10), N2 = c(8, 10),
       alpha = 0.1, Test = 'Strat-Fisher', ss.method = 'exact', N.max = 48)
  }
  a <- ex(0.3)
  b <- ex(0)
  expect_equal(c(a$power.BSSR, a$E.N, b$power.BSSR, b$E.N),
               c(0.567315336728504, 32.6818367739194, 0.0543980782427035, 32.6534657362279),
               tolerance = 1e-10)
})

test_that("StratPowerBSSR agrees with the reference values of the pooled estimate", {
  one <- pb(p = c(0.4, 0.6), Delta.T = 0.3, ss.strata = 'pooled')
  two <- pb(p = c(0.4, 0.6), Delta.T = 0, ss.strata = 'pooled')
  expect_equal(c(one$power.BSSR, one$E.N, two$power.BSSR, two$E.N),
               c(0.811374950528284, 65.21175507318, 0.0517254447539776, 64.9123712443604),
               tolerance = 1e-10)
  ex <- function(dT) {
    pb(p = c(0.4, 0.6), Delta.A = 0.4, Delta.T = dT, N1 = c(8, 10), N2 = c(8, 10),
       alpha = 0.1, Test = 'Unstrat-Fisher', ss.method = 'exact', N.max = 48,
       ss.strata = 'pooled')
  }
  a <- ex(0.3)
  b <- ex(0)
  expect_equal(c(a$power.BSSR, a$E.N, b$power.BSSR, b$E.N),
               c(0.584748568726719, 33.6281430484303, 0.0558427838549891, 33.5446561534546),
               tolerance = 1e-10)
  # With the pooled estimate the final sizes depend only on the total number of interim
  # responders
  expect_identical(attr(one, 'ss.strata'), 'pooled')
  m <- attr(one, 'reestimation')
  expect_true(all(tapply(m$N, m$s1 + m$s2, function(v) length(unique(v))) == 1))
  expect_gt(length(unique(m$N)), 1)
  expect_identical(attr(pb(p = c(0.4, 0.6), Delta.T = 0.3), 'ss.strata'), 'stratified')
})

test_that("StratPowerBSSR returns the fixed design and the distribution of N", {
  res <- pb(p = rbind(c(0.4, 0.6), c(0.1, 0.6), c(0.5, 0.5)), Delta.T = 0.3)
  # p = 0.1 gives p2 = -0.05 in stratum 1 and the scenario is omitted
  expect_identical(nrow(res), 2L)
  expect_equal(res$power.TRAD, StratPower(rbind(c(0.55, 0.75), c(0.65, 0.65)),
                                          rbind(c(0.25, 0.45), c(0.35, 0.35)),
                                          c(10, 10), c(10, 10), 0.05, 'CMH')$Power)
  d <- attr(res, 'N.dist')
  for (k in 1:2) {
    expect_equal(sum(d$prob[d$scenario == k]), 1)
    expect_equal(sum(d$N[d$scenario == k] * d$prob[d$scenario == k]), res$E.N[k])
  }
  expect_s3_class(res, 'bbssr_spowerbssr')
  expect_named(res, c('p.S1', 'p.S2', 'p1.S1', 'p1.S2', 'p2.S1', 'p2.S2', 'power.BSSR',
                      'power.TRAD', 'E.N'))
  expect_identical(attr(res, 'n1.interim'), c(5L, 5L))
  expect_identical(nrow(attr(res, 'reestimation')), 121L)
  expect_error(pb(p = c(0.1, 0.95), Delta.T = 0.3), 'no scenario')
  expect_error(pb(Delta.T = c(0.1, 0.2, 0.3)), 'Delta.T')
})

test_that("StratPowerBSSR simulates designs with three strata", {
  set.seed(5)
  res <- StratPowerBSSR(p = c(0.4, 0.5, 0.6), Delta.A = 0.3, Delta.T = 0.3,
                        N1 = c(8, 8, 8), N2 = c(8, 8, 8), omega = 0.5, r = 1,
                        alpha = 0.05, tar.power = 0.8, Test = 'CMH',
                        ss.method = 'standard', nsim = 2000)
  expect_true('se.BSSR' %in% names(res))
  expect_equal(res$se.BSSR, sqrt(res$power.BSSR * (1 - res$power.BSSR) / 2000))
  expect_identical(attr(res, 'nsim'), 2000L)
  expect_null(attr(res, 'reestimation'))
  expect_error(StratPowerBSSR(p = c(0.4, 0.5, 0.6), Delta.A = 0.3, Delta.T = 0.3,
                              N1 = c(8, 8, 8), N2 = c(8, 8, 8), omega = 0.5, r = 1,
                              alpha = 0.05, tar.power = 0.8, Test = 'CMH'), 'standard')
})

test_that("StratPowerBSSR satisfies the symmetry of the groups and the two-sided union", {
  pw <- function(...) {
    pb(p = rbind(c(0.4, 0.6), c(0.3, 0.5)), ...)
  }
  # Exchanging the groups turns 'less' with the effects -0.3 into 'greater' with 0.3 and
  # leaves the distribution of the pooled interim responders unchanged
  greater <- pw(Delta.A = 0.3, Delta.T = 0.3)
  less <- pw(Delta.A = -0.3, Delta.T = -0.3, alternative = 'less')
  expect_equal(less$power.BSSR, greater$power.BSSR, tolerance = 1e-12)
  expect_equal(less$E.N, greater$E.N, tolerance = 1e-12)
  # The two-sided CMH test at 0.1 rejects where one of the one-sided tests at 0.05 rejects,
  # and its re-estimation uses the same quantile as the one-sided re-estimation at 0.05.
  # The lower tail under the effect 0.3 is the upper tail under -0.3 after the exchange
  two <- pw(Delta.A = 0.3, Delta.T = 0.3, alpha = 0.1, alternative = 'two.sided')
  rev <- pw(Delta.A = 0.3, Delta.T = -0.3)
  expect_gt(rev$power.BSSR[1], 0)
  expect_equal(two$power.BSSR, greater$power.BSSR + rev$power.BSSR, tolerance = 1e-12)
  expect_equal(two$E.N, greater$E.N, tolerance = 1e-12)
})

test_that("StratPowerBSSR with a random allocation agrees with a fixed final total", {
  # With N.min = N.max = 30 every interim outcome leads to 30 patients, so the design is the
  # fixed design of 30 patients with the same allocation, whose power is exact; the average
  # over the sequences must agree with it within its Monte Carlo error
  for (alloc in c('stratum', 'group', 'both')) {
    set.seed(11)
    res <- StratPowerBSSR(p = c(0.4, 0.6), Delta.A = 0.3, Delta.T = 0.3, N1 = c(6, 6),
                          N2 = c(6, 6), omega = 0.5, r = 1, alpha = 0.05, tar.power = 0.8,
                          Test = 'CMH', ss.method = 'standard', N.min = 30, N.max = 30,
                          alloc = alloc, nseq = 40)
    ex <- StratPower(c(0.55, 0.75), c(0.25, 0.45), alpha = 0.05, Test = 'CMH',
                     alloc = alloc, N = 30, frac = c(0.5, 0.5))$Power
    expect_lt(abs(res$power.BSSR - ex), 4 * res$se.BSSR + 1e-12)
    expect_gt(res$se.BSSR, 0)
    expect_equal(res$E.N, 30)
    expect_identical(attr(res, 'N.dist')$N, 30L)
    # The fixed design has the planned total of 24 and the same allocation
    expect_equal(res$power.TRAD,
                 StratPower(c(0.55, 0.75), c(0.25, 0.45), alpha = 0.05, Test = 'CMH',
                            alloc = alloc, N = 24, frac = c(0.5, 0.5))$Power)
  }
})

test_that("StratPowerBSSR with a random allocation is reproducible and keeps settings", {
  run <- function() {
    StratPowerBSSR(p = rbind(c(0.4, 0.6), c(0.5, 0.5)), Delta.A = 0.3, Delta.T = 0.3,
                   N1 = c(6, 6), N2 = c(6, 6), omega = 0.5, r = 1, alpha = 0.05,
                   tar.power = 0.8, Test = 'Strat-Fisher', ss.method = 'standard',
                   N.max = 40, alloc = 'stratum', nseq = 10)
  }
  set.seed(5)
  a <- run()
  set.seed(5)
  b <- run()
  expect_identical(a$power.BSSR, b$power.BSSR)
  expect_named(a, c('p.S1', 'p.S2', 'p1.S1', 'p1.S2', 'p2.S1', 'p2.S2', 'power.BSSR',
                    'power.TRAD', 'E.N', 'se.BSSR'))
  expect_identical(attr(a, 'nseq'), 10L)
  expect_identical(attr(a, 'n.interim'), 12L)
  expect_identical(attr(a, 'block'), 4L)
  expect_null(attr(a, 'reestimation'))
  d <- attr(a, 'N.dist')
  expect_equal(as.vector(tapply(d$prob, d$scenario, sum)), c(1, 1))
  expect_true(all(d$N >= 12 & d$N <= 40))
  expect_equal(a$E.N, as.vector(tapply(d$prob * d$N, d$scenario, sum)))
})

test_that("StratPowerBSSR checks the settings of a random allocation", {
  args <- list(p = c(0.4, 0.6), Delta.A = 0.3, Delta.T = 0.3, N1 = c(6, 6), N2 = c(6, 6),
               omega = 0.5, r = 1, alpha = 0.05, tar.power = 0.8, Test = 'CMH',
               ss.method = 'standard', alloc = 'stratum')
  expect_error(do.call(StratPowerBSSR, c(args, list(nseq = 1))), 'nseq')
  expect_error(do.call(StratPowerBSSR, c(args, list(N.max = 10))), 'N.max')
  expect_error(do.call(StratPowerBSSR, c(args, list(block = 3))), 'whole numbers')
  args3 <- modifyList(args, list(p = c(0.4, 0.5, 0.6), N1 = c(6, 6, 6), N2 = c(6, 6, 6)))
  expect_error(do.call(StratPowerBSSR, args3), 'two strata')
})
