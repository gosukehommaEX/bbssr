tie <- function(...) {
  args <- list(Delta.A = 0.3, N1 = c(10, 10), N2 = c(10, 10), omega = 0.5, r = 1,
               alpha = 0.05, tar.power = 0.8, Test = 'CMH', ss.method = 'standard')
  new <- list(...)
  args[names(new)] <- new
  do.call(StratTypeIErrorBSSR, args)
}

test_that("StratTypeIErrorBSSR agrees with the reference values", {
  res <- tie(theta = rbind(c(0.1, 0.3), c(0.5, 0.5), c(0.9, 0.2)), maximize = 'grid')
  expect_equal(res$TIE.BSSR, c(0.0476476903414493, 0.0512792925484029, 0.0431103731434614),
               tolerance = 1e-10)
  expect_equal(res$TIE.TRAD, c(0.0511555543032792, 0.0445454693253851, 0.0529149719952471),
               tolerance = 1e-10)
  # Largest values over the unit square: grid with local optimization in the reference
  cert <- tie(theta = seq(0, 1, by = 0.1))
  mx <- attr(cert, 'max')
  expect_equal(mx$TIE, c(0.051279292548403, 0.0578880310058594), tolerance = 1e-8)
  expect_true(all(mx$bound >= mx$TIE & mx$bound - mx$TIE <= 1e-12))
  expect_true(all(mx$TIE >= c(max(cert$TIE.BSSR), max(cert$TIE.TRAD)) - 1e-15))
})

test_that("StratTypeIErrorBSSR agrees with the reference values of an unstratified test", {
  # The allocation ratio is 1 in stratum 1 and 3 in stratum 2, so the table of totals
  # mixes strata with different response probabilities in different proportions
  res <- tie(Delta.A = 0.4, N1 = c(6, 9), N2 = c(6, 3), r = c(1, 3), Test = 'Unstrat-Chisq',
             N.max = 72, theta = rbind(c(0.1, 0.3), c(0.5, 0.5), c(0.2, 0.8)),
             maximize = 'grid')
  expect_equal(res$TIE.BSSR, c(0.0598976228961289, 0.0549331938924263, 0.12888549601319),
               tolerance = 1e-10)
  expect_equal(res$TIE.TRAD, c(0.0875167736767266, 0.0505431294441224, 0.126121635879333),
               tolerance = 1e-10)
  expect_identical(attr(res, 'Test'), 'Unstrat-Chisq')
})

test_that("StratTypeIErrorBSSR returns the grid and the settings", {
  res <- tie(theta = c(0.2, 0.5, 0.8))
  expect_s3_class(res, 'bbssr_stie')
  expect_named(res, c('theta.S1', 'theta.S2', 'TIE.BSSR', 'TIE.TRAD'))
  expect_identical(nrow(res), 9L)
  expect_equal(attr(res, 'interval'),
               matrix(c(0.2, 0.2, 0.8, 0.8), 2,
                      dimnames = list(c('theta.S1', 'theta.S2'), NULL)))
  # The fixed design agrees with StratPower under the null hypothesis
  expect_equal(res$TIE.TRAD, StratPower(as.matrix(res[, 1:2]), as.matrix(res[, 1:2]),
                                        c(10, 10), c(10, 10), 0.05, 'CMH')$Power,
               tolerance = 1e-13)
  # The stratified exact test keeps its level in the fixed design
  ex <- tie(Test = 'Strat-Fisher', theta = seq(0, 1, by = 0.25))
  expect_lte(attr(ex, 'max')$bound[2], 0.05)
  expect_error(tie(theta = rbind(c(0.2, 0.3))), 'span an interval')
})

test_that("StratTypeIErrorBSSR simulates designs with three strata", {
  set.seed(3)
  res <- tie(N1 = c(6, 6, 6), N2 = c(6, 6, 6), theta = c(0.3, 0.6), maximize = 'grid',
             nsim = 1000)
  expect_named(res, c('theta.S1', 'theta.S2', 'theta.S3', 'TIE.BSSR', 'TIE.TRAD',
                      'se.BSSR'))
  expect_identical(nrow(res), 2L)
  expect_equal(res$TIE.TRAD, StratPower(rbind(rep(0.3, 3), rep(0.6, 3)),
                                        rbind(rep(0.3, 3), rep(0.6, 3)), c(6, 6, 6),
                                        c(6, 6, 6), 0.05, 'CMH')$Power)
  expect_true(all(is.na(attr(res, 'max')$bound)))
})

test_that("StratTypeIErrorBSSR satisfies the symmetry of the groups and the two-sided rule", {
  th <- rbind(c(0.1, 0.3), c(0.5, 0.5), c(0.9, 0.2))
  g <- tie(theta = th, maximize = 'grid')
  l <- tie(theta = th, maximize = 'grid', Delta.A = -0.3, alternative = 'less')
  expect_equal(l$TIE.BSSR, g$TIE.BSSR, tolerance = 1e-12)
  expect_equal(l$TIE.TRAD, g$TIE.TRAD, tolerance = 1e-12)
  # Under the null hypothesis the two tails have the same probability with equal groups
  two <- tie(theta = th, maximize = 'grid', alpha = 0.1, alternative = 'two.sided')
  expect_equal(two$TIE.BSSR, 2 * g$TIE.BSSR, tolerance = 1e-12)
  expect_equal(two$TIE.TRAD, 2 * g$TIE.TRAD, tolerance = 1e-12)
})

test_that("StratTypeIErrorBSSR simulates three strata in agreement with the exact rate", {
  # With N.min = N.max = 36 every final design is the planned one, so the simulated rate
  # of the re-estimation design estimates the exact rate of the fixed design
  set.seed(4)
  res <- tie(N1 = c(6, 6, 6), N2 = c(6, 6, 6), theta = c(0.3, 0.6), N.min = 36,
             N.max = 36, nsim = 2e4)
  expect_identical(attr(res, 'maximize'), 'grid')
  expect_equal(res$se.BSSR, sqrt(res$TIE.BSSR * (1 - res$TIE.BSSR) / 2e4))
  expect_true(all(abs(res$TIE.BSSR - res$TIE.TRAD) <= 4 * res$se.BSSR))
  expect_error(tie(N1 = c(6, 6, 6), N2 = c(6, 6, 6), maximize = 'certified'),
               'two strata')
})

test_that("StratTypeIErrorBSSR with a random allocation agrees with a fixed final total", {
  # A fixed final total of 30, as in the test of StratPowerBSSR: the average rate over the
  # sequences agrees with the exact rate of the fixed design of 30 patients
  set.seed(12)
  res <- StratTypeIErrorBSSR(Delta.A = 0.3, N1 = c(6, 6), N2 = c(6, 6), omega = 0.5, r = 1,
                             alpha = 0.05, tar.power = 0.8, Test = 'CMH',
                             ss.method = 'standard', N.min = 30, N.max = 30,
                             theta = c(0.3, 0.6), alloc = 'both', nseq = 40)
  pts <- as.matrix(res[, c('theta.S1', 'theta.S2')])
  ex <- StratPower(pts, pts, alpha = 0.05, Test = 'CMH', alloc = 'both', N = 30,
                   frac = c(0.5, 0.5))$Power
  expect_true(all(abs(res$TIE.BSSR - ex) < 4 * res$se.BSSR + 1e-12))
  # The fixed design with the planned total of 24 is exact
  tr <- StratPower(pts, pts, alpha = 0.05, Test = 'CMH', alloc = 'both', N = 24,
                   frac = c(0.5, 0.5))$Power
  expect_equal(res$TIE.TRAD, tr, tolerance = 1e-10)
  m <- attr(res, 'max')
  expect_gte(m$TIE[m$Design == 'BSSR'], max(res$TIE.BSSR))
  expect_gte(m$bound[1], m$TIE[1])
  expect_equal(m$se[1], res$se.BSSR[which.max(res$TIE.BSSR)])
  expect_identical(attr(res, 'nseq'), 40L)
})

test_that("StratTypeIErrorBSSR checks the settings of a random allocation", {
  args <- list(Delta.A = 0.3, N1 = c(6, 6), N2 = c(6, 6), omega = 0.5, r = 1, alpha = 0.05,
               tar.power = 0.8, Test = 'CMH', ss.method = 'standard', alloc = 'stratum')
  expect_error(do.call(StratTypeIErrorBSSR, c(args, list(nseq = 1))), 'nseq')
  args3 <- modifyList(args, list(N1 = c(6, 6, 6), N2 = c(6, 6, 6), alloc = 'group'))
  expect_error(do.call(StratTypeIErrorBSSR, args3), 'two strata')
})

test_that("StratTypeIErrorBSSR simulates a random allocation with three strata", {
  set.seed(22)
  res <- StratTypeIErrorBSSR(Delta.A = 0.3, N1 = c(2, 2, 1), N2 = c(2, 2, 1), omega = 0.5,
                             r = 1, alpha = 0.05, tar.power = 0.8, Test = 'CMH',
                             ss.method = 'standard', frac = c(0.3, 0.3, 0.4), N.max = 20,
                             theta = c(0.3, 0.6), alloc = 'stratum', nsim = 2000)
  expect_named(res, c('theta.S1', 'theta.S2', 'theta.S3', 'TIE.BSSR', 'TIE.TRAD',
                      'se.BSSR', 'se.TRAD'))
  # The same draws as strat_alloc_simulate itself, one point after the other
  a <- check_strat_args(rep(1, 3), rep(1, 3), 0.05, 'CMH', 'greater', 'minlike')
  al <- check_strat_alloc('stratum', NULL, c(0.3, 0.3, 0.4), 1, 4, 3)
  re.args <- list(Delta.A = 0.3, effect = 'RD', frac = c(0.3, 0.3, 0.4), r = c(1, 1, 1),
                  ss.method = 'standard', ss.Test = 'CMH', ss.alpha = 0.05,
                  tar.power = 0.8, search = 'crossing', search.limit = c(2, 50),
                  ss.strata = 'stratified')
  set.seed(22)
  sims <- lapply(c(0.3, 0.6), function(th) {
    strat_alloc_simulate(6, 10, al, 'stratum', rep(th, 3), rep(th, 3), 2000, 0.05, a,
                         re.args, 'observed', FALSE, NULL, 20)
  })
  expect_equal(res$TIE.BSSR, vapply(sims, function(s) mean(s$reject), numeric(1)))
  expect_equal(res$TIE.TRAD, vapply(sims, function(s) mean(s$reject.trad), numeric(1)))
  expect_equal(res$se.TRAD, sqrt(res$TIE.TRAD * (1 - res$TIE.TRAD) / 2000))
  m <- attr(res, 'max')
  expect_equal(m$se, c(res$se.BSSR[which.max(res$TIE.BSSR)],
                       res$se.TRAD[which.max(res$TIE.TRAD)]))
  expect_identical(attr(res, 'maximize'), 'grid')
  expect_null(attr(res, 'nseq'))
})
