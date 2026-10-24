cr_map <- function() {
  strat_map(c(3L, 3L), c(3L, 3L), 0.05, 'CMH', 'greater', 'minlike',
            list(Delta.A = 0.3, effect = 'RD', frac = c(0.5, 0.5), r = c(1, 1),
                 ss.method = 'standard', ss.Test = 'CMH', ss.alpha = 0.05,
                 tar.power = 0.8, N1.plan = c(6L, 6L), N2.plan = c(6L, 6L),
                 restricted = FALSE, N.min = NULL, N.max = NULL, search = 'crossing',
                 search.limit = c(2, 50)))
}

test_that("strat_cond_reject returns the summaries per interim outcome", {
  map <- cr_map()
  cr <- strat_cond_reject(c(3L, 3L), c(3L, 3L), map$designs, map$id, 0.05, FALSE,
                          c(0.2, 0.6), c(0.4, 0.6))
  expect_length(cr$crp.max, 49)
  expect_identical(dim(cr$tie.s), c(49L, 2L))
  expect_length(cr$crp, 0)
  expect_true(all(cr$crp.min <= cr$crp.max))
  expect_true(all(cr$exceed >= 0))
  # Averaged over the interim outcomes, the rates given s give the type I error rate
  prob <- outer(stats::dbinom(0:6, 6, 0.2), stats::dbinom(0:6, 6, 0.4))
  expect_equal(sum(as.vector(prob) * cr$tie.s[, 1]),
               strat_bssr_power(c(3L, 3L), c(3L, 3L), map$designs, map$id,
                                rbind(c(0.2, 0.4)), rbind(c(0.2, 0.4))),
               tolerance = 1e-13)
  none <- strat_cond_reject(c(3L, 3L), c(3L, 3L), map$designs, map$id, 0.05, TRUE,
                            numeric(0), numeric(0))
  expect_identical(dim(none$tie.s), c(49L, 0L))
  expect_gt(length(none$crp), 49)
  expect_identical(length(none$crp), length(none$w2))
})

test_that("strat_cond_reject checks its arguments", {
  map <- cr_map()
  f <- function(alpha = 0.05, t1 = 0.5, t2 = 0.5) {
    strat_cond_reject(c(3L, 3L), c(3L, 3L), map$designs, map$id, alpha, FALSE, t1, t2)
  }
  expect_error(f(t2 = numeric(0)), 'same length')
  expect_error(f(t1 = 1.5), 'unit square')
  expect_error(f(alpha = 0), 'alpha')
})
