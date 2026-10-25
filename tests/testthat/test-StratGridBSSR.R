# Reference values from tools/reference/reference_strat.py (blocks_summary_sbssr), an
# independent Python implementation that re-estimates the sample size by the same rules
# and enumerates the interim and the second-stage outcomes directly
test_that("StratGridBSSR agrees with an independent implementation", {
  design <- data.frame(Test = c('CMH', 'Strat-Fisher'), ss.method = c('standard', 'exact'),
                       Delta.A = c(0.3, 0.4), alpha = c(0.05, 0.1), N.max = c(NA, 48),
                       N1.S1 = c(10, 8), N1.S2 = c(10, 10), N2.S1 = c(10, 8),
                       N2.S2 = c(10, 10))
  res <- StratGridBSSR(design, p = c(0.4, 0.6), Delta.T = 0.3, omega = 0.5, r = 1,
                       tar.power = 0.8)
  expect_s3_class(res, 'bbssr_sgrid')
  expect_s3_class(res, 'bbssr_grid')
  cols <- c('power.BSSR', 'E.N', 'SD.N', 'N.q25', 'N.q50', 'N.q75', 'P.N.max')
  got <- c(unlist(res[1, cols]), unlist(res[2, cols]))
  expect_equal(unname(got),
               c(0.770441457550367, 59.740212490731, 7.12815613980113, 56, 60, 64,
                 0.198521899360162, 0.567315336728504, 32.6818367739194, 2.4998983667878,
                 32, 34, 34, 5.90356207233369e-05),
               tolerance = 1e-10)
  expect_equal(res$design, c(1, 2))
  expect_equal(res$Test, c('CMH', 'Strat-Fisher'))
  expect_equal(res$N1.S1, c(10, 8))
  expect_equal(res$n1.interim.S1, c(5, 4))
  expect_equal(res$n2.interim.S2, c(5, 5))
  expect_equal(res$Delta.T.S2, c(0.3, 0.3))
  expect_identical(res$scenario, c(1L, 1L))
  expect_equal(res$p1.S1, c(0.55, 0.55))
  expect_null(res$se.BSSR)
  expect_true(is.na(res$N.max[1]))
})

test_that("each design gives the result of StratPowerBSSR and its summary", {
  # The second design gives the interim sizes instead of omega. The second scenario has
  # p2 = -0.05 in stratum 1 and is dropped
  design <- data.frame(Test = c('CMH', 'Cochran'), omega = c(0.5, NA),
                       n1.interim.S1 = c(NA, 4), n1.interim.S2 = c(NA, 6),
                       n2.interim.S1 = c(NA, 4), n2.interim.S2 = c(NA, 6))
  p <- rbind(c(0.4, 0.6), c(0.1, 0.6), c(0.5, 0.5))
  res <- StratGridBSSR(design, p = p, Delta.A = 0.3, N1 = c(10, 10), N2 = c(10, 10), r = 1,
                       alpha = 0.05, tar.power = 0.8, ss.method = 'standard')
  direct <- list(
    StratPowerBSSR(p = p, Delta.A = 0.3, Delta.T = 0.3, N1 = c(10, 10), N2 = c(10, 10),
                   omega = 0.5, r = 1, alpha = 0.05, tar.power = 0.8, Test = 'CMH',
                   ss.method = 'standard'),
    StratPowerBSSR(p = p, Delta.A = 0.3, Delta.T = 0.3, N1 = c(10, 10), N2 = c(10, 10),
                   n.interim = rbind(c(4, 6), c(4, 6)), r = 1, alpha = 0.05,
                   tar.power = 0.8, Test = 'Cochran', ss.method = 'standard')
  )
  for (i in 1:2) {
    k <- res$design == i
    d <- direct[[i]]
    s <- summary(d)
    expect_identical(res$scenario[k], c(1L, 3L))
    expect_identical(res$p.S2[k], d$p.S2)
    expect_identical(res$p1.S1[k], d$p1.S1)
    expect_identical(res$power.BSSR[k], d$power.BSSR)
    expect_identical(res$power.TRAD[k], d$power.TRAD)
    expect_identical(res$E.N[k], d$E.N)
    expect_identical(res$SD.N[k], s$SD.N)
    expect_identical(res$N.q75[k], s$N.q75)
    expect_identical(res$P.N.max[k], s$P.N.max)
    expect_equal(unique(res$n1.interim.S2[k]), attr(d, 'n1.interim')[2])
  }
  d <- attr(res, 'designs')
  expect_equal(nrow(d), 2)
  expect_named(d, c('design', 'Test', 'omega', 'N1.S1', 'N1.S2', 'N2.S1', 'N2.S2',
                    'n1.interim.S1', 'n1.interim.S2', 'n2.interim.S1', 'n2.interim.S2',
                    'Delta.T.S1', 'Delta.T.S2'))
  # summary.bbssr_grid applies to the stratified grid
  sm <- summary(res)
  expect_equal(sm$power.BSSR.max, c(max(direct[[1]]$power.BSSR),
                                    max(direct[[2]]$power.BSSR)))
})

test_that("the planned group sizes follow from p.plan", {
  # Pooled planning probabilities 0.3 and 0.5, risk difference 0.4, fractions (0.4, 0.6)
  design <- data.frame(p.plan.S1 = 0.3, p.plan.S2 = 0.5)
  res <- StratGridBSSR(design, p = c(0.3, 0.5), Delta.A = 0.4, frac = c(0.4, 0.6),
                       omega = 0.5, r = 1, alpha = 0.025, tar.power = 0.8, Test = 'CMH',
                       ss.method = 'standard')
  d <- attr(res, 'designs')
  expect_equal(c(d$N1.S1, d$N1.S2, d$N2.S1, d$N2.S2),
               c(9, 14, 9, 14))
  ss <- StratSampleSize(c(0.5, 0.7), c(0.1, 0.3), frac = c(0.4, 0.6), r = 1,
                        alpha = 0.025, tar.power = 0.8, Test = 'CMH', method = 'standard')
  expect_equal(c(d$N1.S1, d$N1.S2), ss$N1)
  expect_equal(c(d$N2.S1, d$N2.S2), ss$N2)
  expect_equal(d$p.plan.S2, 0.5)
  expect_error(StratGridBSSR(design, p = c(0.3, 0.5), Delta.A = 0.4, omega = 0.5, r = 1,
                             alpha = 0.025, tar.power = 0.8, Test = 'CMH',
                             ss.method = 'standard'),
               'design 1: frac must be given with p.plan')
})

test_that("arguments of the strata are taken from the columns and the shared vectors", {
  run <- function(design, ...) {
    StratGridBSSR(design, p = c(0.4, 0.6), Delta.A = 0.3, N1 = c(10, 10), N2 = c(10, 10),
                  omega = 0.5, alpha = 0.05, tar.power = 0.8, Test = 'CMH',
                  ss.method = 'standard', ...)
  }
  a <- run(data.frame(r.S1 = 1, r.S2 = 1, Delta.T.S1 = 0.3, Delta.T.S2 = 0))
  b <- run(data.frame(r = 1), Delta.T = c(0.3, 0))
  c0 <- StratPowerBSSR(p = c(0.4, 0.6), Delta.A = 0.3, Delta.T = c(0.3, 0),
                       N1 = c(10, 10), N2 = c(10, 10), omega = 0.5, r = 1, alpha = 0.05,
                       tar.power = 0.8, Test = 'CMH', ss.method = 'standard')
  expect_identical(a$power.BSSR, c0$power.BSSR)
  expect_identical(b$power.BSSR, c0$power.BSSR)
  expect_equal(c(a$Delta.T.S1, a$Delta.T.S2), c(0.3, 0))
  expect_equal(a$r.S2, 1)
  # The assumed effect is the default true effect
  expect_equal(run(data.frame(r = 1))$Delta.T.S2, 0.3)
  expect_error(run(data.frame(r = 1, r.S1 = 1, r.S2 = 1)), 'with and without the suffix')
  expect_error(run(data.frame(r.S1 = 1)), 'number the strata 1 to 2')
  expect_error(run(data.frame(r.S1 = c(1, 1), r.S2 = c(1, NA))),
               'design 2: the columns r.S1, ... are missing for some strata')
  expect_error(run(data.frame(r = 1, foo = 1)), 'unknown argument')
  expect_error(run(data.frame(r = 1, Test = 'CMH')), 'given both in design and as a shared')
  expect_error(run(data.frame(r = 1, Delta.T = 0.3), Delta.T = 0.3), 'give Delta.T either')
  expect_error(run(data.frame(r = 1), Delta.T = c(0.3, 0.2, 0.1)),
               'design 1: Delta.T must be a single value or one value per stratum')
  expect_error(StratGridBSSR(data.frame(r = 1), p = 0.4), 'two or more strata')
  expect_error(StratGridBSSR(data.frame(r = 1), p = c(0.4, NA)), 'without missing values')
  expect_error(StratGridBSSR(data.frame(), p = c(0.4, 0.6)), 'at least one row')
  expect_error(StratGridBSSR(data.frame(r = 1), p = c(0.4, 0.6), 0.3, FALSE, FALSE, 1),
               'must be named')
})

test_that("type1 adds the largest type I error rates of StratTypeIErrorBSSR", {
  args <- list(Delta.A = 0.3, N1 = c(6, 6), N2 = c(6, 6), omega = 0.5, r = 1,
               alpha = 0.05, tar.power = 0.8, Test = 'CMH', ss.method = 'standard',
               N.max = 40)
  res <- do.call(StratGridBSSR, c(list(design = data.frame(Test = 'CMH'),
                                       p = c(0.4, 0.6), type1 = TRUE), args[-8]))
  m <- attr(do.call(StratTypeIErrorBSSR, args), 'max')
  d <- attr(res, 'designs')
  expect_identical(c(d$TIE.BSSR, d$TIE.TRAD), m$TIE)
  expect_identical(c(d$theta.S1.BSSR, d$theta.S2.BSSR), unlist(m[1, c('theta.S1',
                                                                      'theta.S2')],
                                                                use.names = FALSE))
  expect_identical(c(d$theta.S1.TRAD, d$theta.S2.TRAD), unlist(m[2, c('theta.S1',
                                                                      'theta.S2')],
                                                                use.names = FALSE))
  expect_identical(c(d$bound.BSSR, d$bound.TRAD), m$bound)
  expect_null(d$se.BSSR)
  expect_true(attr(res, 'type1'))
})

test_that("StratGridBSSR names the design of an error or a warning", {
  design <- data.frame(Test = c('CMH', 'foo'))
  expect_error(StratGridBSSR(design, p = c(0.4, 0.6), Delta.A = 0.3, N1 = c(10, 10),
                             N2 = c(10, 10), omega = 0.5, r = 1, alpha = 0.05,
                             tar.power = 0.8, ss.method = 'standard'),
               'design 2: ')
  # N1 = 11 is not ceiling(r N2) in stratum 1
  design <- data.frame(N1.S1 = 11, N1.S2 = 10, N2.S1 = 10, N2.S2 = 10)
  expect_warning(
    res <- StratGridBSSR(design, p = c(0.4, 0.6), Delta.A = 0.3, omega = 0.5, r = 1,
                         alpha = 0.05, tar.power = 0.8, Test = 'CMH',
                         ss.method = 'standard'),
    'design 1: N1 differs from ceiling')
  # A design with no column other than the group sizes
  expect_named(attr(res, 'designs'), c('design', 'N1.S1', 'N1.S2', 'N2.S1', 'N2.S2',
                                       'n1.interim.S1', 'n1.interim.S2', 'n2.interim.S1',
                                       'n2.interim.S2', 'Delta.T.S1', 'Delta.T.S2'))
  expect_message(StratGridBSSR(data.frame(Test = 'Cochran'), p = c(0.4, 0.6),
                               Delta.A = 0.3, N1 = c(10, 10), N2 = c(10, 10), omega = 0.5,
                               r = 1, alpha = 0.05, tar.power = 0.8,
                               ss.method = 'standard', verbose = TRUE),
                 'design 1 of 1')
})

test_that("StratGridBSSR evaluates random allocations reproducibly", {
  run <- function() {
    set.seed(8)
    StratGridBSSR(data.frame(alloc = c('fixed', 'stratum')), p = c(0.4, 0.6),
                  Delta.A = 0.3, N1 = c(6, 6), N2 = c(6, 6), omega = 0.5, r = 1,
                  alpha = 0.05, tar.power = 0.8, Test = 'CMH', ss.method = 'standard',
                  nseq = 5, N.max = 40)
  }
  res <- run()
  expect_identical(run(), res)
  expect_true(is.na(res$se.BSSR[1]))
  expect_false(is.na(res$se.BSSR[2]))
  expect_null(res$se.TRAD)
  # The first design draws no random numbers
  set.seed(8)
  d <- StratPowerBSSR(p = c(0.4, 0.6), Delta.A = 0.3, Delta.T = 0.3, N1 = c(6, 6),
                      N2 = c(6, 6), omega = 0.5, r = 1, alpha = 0.05, tar.power = 0.8,
                      Test = 'CMH', ss.method = 'standard', N.max = 40, alloc = 'stratum',
                      nseq = 5)
  expect_identical(res$power.BSSR[2], d$power.BSSR)
  expect_identical(res$se.BSSR[2], d$se.BSSR)
  expect_identical(res$P.N.max[2], summary(d)$P.N.max)
})
