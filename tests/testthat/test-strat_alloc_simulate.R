re_args <- function(frac) {
  list(Delta.A = 0.3, effect = 'RD', frac = frac, r = rep(1, length(frac)),
       ss.method = 'standard', ss.Test = 'CMH', ss.alpha = 0.05, tar.power = 0.8,
       search = 'crossing', search.limit = c(2, 50), ss.strata = 'stratified')
}
a_cmh <- function(J) {
  check_strat_args(rep(1, J), rep(1, J), 0.05, 'CMH', 'greater', 'minlike')
}

test_that("strat_alloc_simulate agrees with the exact rates with three strata", {
  # The final analysis has the first 10 patients, whatever the interim outcome
  # (N.min = N.max = 10), and the interim analysis the first 6, drawn again until every
  # stratum has one of them; the fixed design has 10 patients. Exact values from
  # tools/reference/reference_strat.py (alloc_power_strata), for the true probabilities
  # (0.55, 0.65, 0.75) and (0.25, 0.35, 0.45), and theta = 0.3 and 0.6 in every stratum
  frac <- c(0.3, 0.3, 0.4)
  sc <- list(list(c(0.55, 0.65, 0.75), c(0.25, 0.35, 0.45)), list(rep(0.3, 3), rep(0.3, 3)),
             list(rep(0.6, 3), rep(0.6, 3)))
  exact <- list(
    stratum = list(bssr = c(0.230776654915696, 0.0449512053534727, 0.0510778201204156),
                   trad = c(0.230262063558603, 0.0444657581797135, 0.050859908031419)),
    both = list(bssr = c(0.190561177909111, 0.0408275669101311, 0.0450138837018902),
                trad = c(0.190441306733802, 0.0403084968958065, 0.0445455664549969)))
  nsim <- 20000
  for (alloc in c('stratum', 'both')) {
    al <- check_strat_alloc(alloc, NULL, frac, 1, 4, 3)
    set.seed(31)
    for (k in 1:3) {
      sim <- strat_alloc_simulate(6, 10, al, alloc, sc[[k]][[1]], sc[[k]][[2]], nsim, 0.05,
                                  a_cmh(3), re_args(frac), 'observed', FALSE, 10, 10)
      expect_true(all(sim$N == 10))
      for (d in c('bssr', 'trad')) {
        ex <- exact[[alloc]][[d]][k]
        got <- mean(if (d == 'bssr') sim$reject else sim$reject.trad)
        expect_lt(abs(got - ex), 4 * sqrt(ex * (1 - ex) / nsim),
                  label = paste(alloc, d, k))
      }
    }
  }
})

test_that("strat_alloc_simulate with two strata agrees with the exact computations", {
  frac <- c(0.5, 0.5)
  p1 <- c(0.55, 0.75)
  p2 <- c(0.25, 0.45)
  nsim <- 20000
  se <- function(p) sqrt(p * (1 - p) / nsim)
  for (alloc in c('stratum', 'both')) {
    al <- check_strat_alloc(alloc, NULL, frac, 1, 4)
    set.seed(32)
    # A final total fixed at 30 after 16 interim patients, who leave a stratum empty with
    # the probability 2 / 2^16, and the fixed design of 24 patients, against StratPower
    sim <- strat_alloc_simulate(16, 24, al, alloc, p1, p2, nsim, 0.05, a_cmh(2),
                                re_args(frac), 'observed', FALSE, 30, 30)
    ex30 <- StratPower(p1, p2, alpha = 0.05, Test = 'CMH', alloc = alloc, N = 30,
                       frac = frac)$Power
    ex24 <- StratPower(p1, p2, alpha = 0.05, Test = 'CMH', alloc = alloc, N = 24,
                       frac = frac)$Power
    expect_lt(abs(mean(sim$reject) - ex30), 4 * se(ex30), label = paste(alloc, 'N 30'))
    expect_lt(abs(mean(sim$reject.trad) - ex24), 4 * se(ex24),
              label = paste(alloc, 'fixed'))
    # A re-estimated total, against the average over the allocation sequences
    sim <- strat_alloc_simulate(12, 24, al, alloc, p1, p2, nsim, 0.05, a_cmh(2),
                                re_args(frac), 'observed', FALSE, NULL, 40)
    expect_true(all(sim$N >= 12 & sim$N <= 40) && length(unique(sim$N)) > 1)
    sq <- StratPowerBSSR(p = c(0.4, 0.6), Delta.A = 0.3, Delta.T = 0.3, N1 = c(6, 6),
                         N2 = c(6, 6), omega = 0.5, r = 1, alpha = 0.05, tar.power = 0.8,
                         Test = 'CMH', ss.method = 'standard', N.max = 40, alloc = alloc,
                         nseq = 100)
    got <- mean(sim$reject)
    expect_lt(abs(got - sq$power.BSSR), 4 * sqrt(se(got)^2 + sq$se.BSSR^2),
              label = paste(alloc, 're-estimated'))
  }
})

test_that("strat_alloc_simulate stops when no interim analysis has every stratum", {
  al <- check_strat_alloc('both', NULL, c(1, 1, 1), 1, 4, 3)
  set.seed(33)
  expect_error(strat_alloc_simulate(2, 10, al, 'both', rep(0.5, 3), rep(0.5, 3), 10, 0.05,
                                    a_cmh(3), re_args(rep(1, 3) / 3), 'observed', FALSE,
                                    NULL, NULL),
               'no patient in one of the strata')
})

test_that("strat_alloc_simulate is reproducible", {
  al <- check_strat_alloc('stratum', NULL, c(0.2, 0.3, 0.5), 1, 4, 3)
  run <- function() {
    strat_alloc_simulate(12, 18, al, 'stratum', c(0.5, 0.6, 0.7), c(0.3, 0.4, 0.5), 200,
                         0.05, a_cmh(3), re_args(c(0.2, 0.3, 0.5)), 'observed', FALSE,
                         NULL, 36)
  }
  set.seed(34)
  x <- run()
  set.seed(34)
  expect_identical(run(), x)
  expect_length(x$reject, 200)
  expect_type(x$reject.trad, 'logical')
  expect_true(all(x$N >= 12 & x$N <= 36))
})
