cells <- rbind(c(0, 0), c(3, 7), c(6, 4), c(8, 10))
flat <- function(re) as.vector(t(cbind(re$N1, re$N2)))

test_that("strat_reestimate agrees with the reference final sizes", {
  std <- strat_reestimate(cells, c(5L, 5L), c(5L, 5L), 0.3, 'RD', c(0.5, 0.5), c(1, 1),
                          'standard', 'CMH', 0.05, 0.8, 'greater', 'minlike', c(10L, 10L),
                          c(10L, 10L), FALSE, NULL, NULL, 'crossing', c(2, 50))
  expect_equal(flat(std), c(5, 5, 5, 5, 14, 14, 14, 14, 16, 16, 16, 16, 5, 5, 5, 5))
  ex <- strat_reestimate(cells, c(4L, 5L), c(4L, 5L), 0.4, 'RD', c(8, 10) / 18, c(1, 1),
                         'exact', 'Strat-Fisher', 0.1, 0.8, 'greater', 'minlike',
                         c(8L, 10L), c(8L, 10L), FALSE, NULL, 48, 'crossing', c(2, 50))
  expect_equal(flat(ex), c(11, 13, 11, 13, 7, 9, 7, 9, 7, 9, 7, 9, 11, 13, 11, 13))
  expect_true(all(is.na(ex$n.raw)))
})

test_that("strat_reestimate applies the bounds and keeps the planned sizes without effect", {
  args <- list(s = cells, n1 = c(5L, 5L), n2 = c(5L, 5L), Delta.A = 0.3, effect = 'RD',
               frac = c(0.5, 0.5), r = c(1, 1), ss.method = 'standard', ss.Test = 'CMH',
               ss.alpha = 0.05, tar.power = 0.8, alternative = 'greater',
               tsmethod = 'minlike', N1.plan = c(10L, 10L), N2.plan = c(10L, 10L),
               restricted = TRUE, N.min = NULL, N.max = 60, search = 'crossing',
               search.limit = c(2, 50))
  re <- do.call(strat_reestimate, args)
  tot <- rowSums(re$N1) + rowSums(re$N2)
  expect_true(all(re$N1 >= 10 & re$N2 >= 10))
  expect_true(all(tot <= 60))
  expect_equal(re$hat.p[2, ], c(0.3, 0.7))
  expect_equal(re$hat.p1[2, ], c(0.45, 0.85))
  # With the odds ratio, a stratum without responders has equal planning probabilities,
  # and when this holds in every stratum the planned sizes are kept
  args$effect <- 'OR'
  args$Delta.A <- 3
  args$restricted <- FALSE
  args$N.max <- NULL
  args$s <- rbind(c(0, 10), c(10, 0), c(0, 0))
  re <- do.call(strat_reestimate, args)
  expect_true(all(is.nan(re$n.raw)))
  expect_identical(re$N1, matrix(10L, 3, 2))
  args$N1.plan <- NULL
  args$N2.plan <- NULL
  expect_error(do.call(strat_reestimate, args), 'supply the planned sizes')
})
