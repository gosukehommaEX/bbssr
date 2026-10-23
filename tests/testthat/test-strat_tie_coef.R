test_that("strat_tie_coef gives the rejection probability under the null hypothesis", {
  map <- strat_map(c(3L, 3L), c(3L, 3L), 0.05, 'CMH', 'greater', 'minlike',
                   list(Delta.A = 0.3, effect = 'RD', frac = c(0.5, 0.5), r = c(1, 1),
                        ss.method = 'standard', ss.Test = 'CMH', ss.alpha = 0.05,
                        tar.power = 0.8, N1.plan = c(6L, 6L), N2.plan = c(6L, 6L),
                        restricted = FALSE, N.min = NULL, N.max = NULL,
                        search = 'crossing', search.limit = c(2, 50)))
  coef <- strat_tie_coef(c(3L, 3L), c(3L, 3L), map$designs, map$id)
  expect_identical(dim(coef), c(max(map$map$N1.S1 + map$map$N2.S1) + 1L,
                                max(map$map$N1.S2 + map$map$N2.S2) + 1L))
  expect_true(all(coef >= -1e-15))
  th <- rbind(c(0.1, 0.8), c(0.5, 0.5), c(0.95, 0.3))
  expect_equal(strat_bern_eval(coef, th[, 1], th[, 2]),
               strat_bssr_power(c(3L, 3L), c(3L, 3L), map$designs, map$id, th, th),
               tolerance = 1e-13)
})
