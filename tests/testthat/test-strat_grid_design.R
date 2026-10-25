base_args <- function(...) {
  a <- list(Delta.A = 0.3, r = 1, alpha = 0.05, tar.power = 0.8, Test = 'CMH',
            ss.method = 'standard', N1 = c(10, 10), N2 = c(10, 10), omega = 0.5)
  new <- list(...)
  a[names(new)] <- new
  a
}

test_that("strat_grid_design evaluates one design with StratPowerBSSR", {
  res <- strat_grid_design(base_args(), c(0.4, 0.6), NULL, FALSE)
  expect_named(res, c('N1', 'N2', 'Delta.T', 'power', 'summary'))
  d <- StratPowerBSSR(p = c(0.4, 0.6), Delta.A = 0.3, Delta.T = 0.3, N1 = c(10, 10),
                      N2 = c(10, 10), omega = 0.5, r = 1, alpha = 0.05, tar.power = 0.8,
                      Test = 'CMH', ss.method = 'standard')
  expect_identical(res$power$power.BSSR, d$power.BSSR)
  expect_identical(res$summary, summary(d))
  # The true effect: the assumed one by default, the one in the list, or the argument,
  # repeated over the strata
  expect_equal(res$Delta.T, c(0.3, 0.3))
  expect_equal(strat_grid_design(base_args(Delta.T = c(0.2, 0.1)), c(0.4, 0.6), NULL,
                                 FALSE)$Delta.T, c(0.2, 0.1))
  expect_equal(strat_grid_design(base_args(Delta.T = 0.2), c(0.4, 0.6), 0.1,
                                 FALSE)$Delta.T, c(0.1, 0.1))
  # Interim sizes given by the two groups
  b <- strat_grid_design(base_args(omega = NULL, n1.interim = c(4, 6),
                                   n2.interim = c(4, 6)), c(0.4, 0.6), NULL, FALSE)
  expect_equal(attr(b$power, 'n1.interim'), c(4, 6))
})

test_that("strat_grid_design checks the arguments of the design", {
  run <- function(a, Delta.T = NULL) strat_grid_design(a, c(0.4, 0.6), Delta.T, FALSE)
  for (nm in c('Delta.A', 'r', 'alpha', 'tar.power', 'Test')) {
    a <- base_args()
    a[[nm]] <- NULL
    expect_error(run(a), paste(nm, 'must be given'))
  }
  expect_error(run(base_args(omega = NULL, n1.interim = c(4, 6))),
               'give both n1.interim and n2.interim')
  a <- base_args()
  a$N2 <- NULL
  expect_error(run(a), 'give both N1 and N2')
  expect_error(run(base_args(p.plan = c(0.3, 0.5))), 'either N1 and N2 or p.plan, not both')
  a$N1 <- NULL
  expect_error(run(a), 'give either N1 and N2 or p.plan')
  expect_error(run(c(a, list(p.plan = c(0.3, 0.5)))), 'frac must be given')
  expect_error(run(c(a, list(p.plan = c(0.3, 0.5), frac = c(1, 1), effect = 'foo'))),
               'effect must be')
  expect_error(run(base_args(), Delta.T = c(0.1, 0.2, 0.3)), 'Delta.T must be')
  expect_error(run(base_args(), Delta.T = NA_real_), 'Delta.T must be')
})

test_that("strat_grid_design plans the group sizes and the type I error rates", {
  a <- base_args(p.plan = c(0.3, 0.5), frac = c(0.4, 0.6), Delta.A = 0.4, alpha = 0.025)
  a$N1 <- NULL
  a$N2 <- NULL
  res <- strat_grid_design(a, c(0.3, 0.5), NULL, FALSE)
  ss <- StratSampleSize(c(0.5, 0.7), c(0.1, 0.3), frac = c(0.4, 0.6), r = 1,
                        alpha = 0.025, tar.power = 0.8, Test = 'CMH', method = 'standard')
  expect_identical(res$N1, ss$N1)
  expect_identical(res$N2, ss$N2)
  # The odds ratio divides the planning probabilities
  o <- strat_grid_design(modifyList(a, list(effect = 'OR', Delta.A = 6)), c(0.3, 0.5),
                         NULL, FALSE)
  sp <- lapply(c(0.3, 0.5), split_pooled, Delta = 6, r = 1, effect = 'OR')
  so <- StratSampleSize(vapply(sp, function(s) s$p1, numeric(1)),
                        vapply(sp, function(s) s$p2, numeric(1)), frac = c(0.4, 0.6),
                        r = 1, alpha = 0.025, tar.power = 0.8, Test = 'CMH',
                        method = 'standard')
  expect_identical(o$N1, so$N1)
  t1 <- strat_grid_design(base_args(N1 = c(6, 6), N2 = c(6, 6), N.max = 40), c(0.4, 0.6),
                          NULL, TRUE)$type1
  expect_named(t1, c('TIE.BSSR', 'TIE.TRAD', 'theta.S1.BSSR', 'theta.S1.TRAD',
                     'theta.S2.BSSR', 'theta.S2.TRAD', 'bound.BSSR', 'bound.TRAD'))
  expect_true(all(t1$bound.BSSR >= t1$TIE.BSSR))
})
