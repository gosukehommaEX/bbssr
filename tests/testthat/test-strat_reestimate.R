cells <- rbind(c(0, 0), c(3, 7), c(6, 4), c(8, 10))
flat <- function(re) as.vector(t(cbind(re$N1, re$N2)))

test_that("strat_reestimate agrees with the reference final sizes", {
  std <- strat_reestimate(cells, c(5L, 5L), c(5L, 5L), 0.3, 'RD', c(0.5, 0.5), c(1, 1),
                          'standard', 'CMH', 0.05, 0.8, 'greater', 'minlike', c(10L, 10L),
                          c(10L, 10L), FALSE, NULL, NULL, 'crossing', c(2, 50))
  expect_equal(flat(std), c(5, 5, 5, 5, 14, 14, 14, 14, 16, 16, 16, 16, 6, 6, 6, 6))
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
  # Assigning list(NULL) keeps the element, which args$N.max <- NULL would remove
  args['N.max'] <- list(NULL)
  args$s <- rbind(c(0, 10), c(10, 0), c(0, 0))
  re <- do.call(strat_reestimate, args)
  expect_true(all(is.nan(re$n.raw)))
  expect_identical(re$N1, matrix(10L, 3, 2))
  args['N1.plan'] <- list(NULL)
  args['N2.plan'] <- list(NULL)
  expect_error(do.call(strat_reestimate, args), 'supply the planned sizes')
})

test_that("strat_reestimate requires a numeric target power", {
  expect_error(strat_reestimate(cells, c(5L, 5L), c(5L, 5L), 0.3, 'RD', c(0.5, 0.5),
                                c(1, 1), 'standard', 'CMH', 0.05, '0.8', 'greater',
                                'minlike', c(10L, 10L), c(10L, 10L), FALSE, NULL, NULL,
                                'crossing', c(2, 50)), 'tar.power')
})

test_that("strat_reestimate agrees with the reference sizes for the risk and odds ratio", {
  re <- function(Delta.A, effect) {
    strat_reestimate(cells, c(5L, 5L), c(5L, 5L), Delta.A, effect, c(0.5, 0.5), c(1, 1),
                     'standard', 'CMH', 0.05, 0.8, 'greater', 'minlike', c(10L, 10L),
                     c(10L, 10L), FALSE, NULL, NULL, 'crossing', c(2, 50))
  }
  # Without interim responders neither ratio shows an effect, and the planned sizes are
  # kept
  expect_equal(flat(re(1.5, 'RR')),
               c(10, 10, 10, 10, 32, 32, 32, 32, 37, 37, 37, 37, 5, 5, 5, 5))
  expect_equal(flat(re(3, 'OR')),
               c(10, 10, 10, 10, 26, 26, 26, 26, 23, 23, 23, 23, 70, 70, 70, 70))
})

test_that("strat_reestimate agrees with the reference sizes of the pooled estimate", {
  re <- function(Test, ss.strata) {
    strat_reestimate(cells, c(5L, 5L), c(5L, 5L), 0.3, 'RD', c(0.5, 0.5), c(1, 1),
                     'standard', Test, 0.05, 0.8, 'greater', 'minlike', c(10L, 10L),
                     c(10L, 10L), FALSE, NULL, NULL, 'crossing', c(2, 50), ss.strata)
  }
  expect_equal(flat(re('CMH', 'pooled')),
               c(5, 5, 5, 5, 17, 17, 17, 17, 17, 17, 17, 17, 7, 7, 7, 7))
  expect_equal(flat(re('Unstrat-Chisq', 'stratified')),
               c(5, 5, 5, 5, 17, 17, 17, 17, 17, 17, 17, 17, 7, 7, 7, 7))
  expect_equal(flat(re('Unstrat-Chisq', 'pooled')),
               c(5, 5, 5, 5, 17, 17, 17, 17, 17, 17, 17, 17, 7, 7, 7, 7))
  # The pooled estimate is the same in every stratum and depends only on the total
  p <- re('CMH', 'pooled')
  expect_equal(p$hat.p[, 1], p$hat.p[, 2])
  expect_equal(p$hat.p[, 1], rowSums(cells) / 20)
  q <- strat_reestimate(rbind(c(3, 7), c(7, 3), c(10, 0)), c(5L, 5L), c(5L, 5L), 0.3, 'RD',
                        c(0.5, 0.5), c(1, 1), 'standard', 'CMH', 0.05, 0.8, 'greater',
                        'minlike', c(10L, 10L), c(10L, 10L), FALSE, NULL, NULL, 'crossing',
                        c(2, 50), 'pooled')
  expect_identical(q$N1[1, ], q$N1[2, ])
  expect_identical(q$N1[1, ], q$N1[3, ])
  expect_error(re('CMH', 'by stratum'))
})

test_that("the exact search of strat_reestimate does not go beyond N.max", {
  # With the assumed odds ratio 6 the normal approximation of some interim outcomes
  # lies within the bound and that of others beyond it
  s <- as.matrix(expand.grid(0:8, 0:8))
  args <- list(s = s, n1 = c(4L, 4L), n2 = c(4L, 4L), Delta.A = 6, effect = 'OR',
               frac = c(0.5, 0.5), r = c(1, 1), ss.method = 'exact', ss.Test = 'CMH',
               ss.alpha = 0.05, tar.power = 0.8, alternative = 'greater',
               tsmethod = 'minlike', N1.plan = c(8L, 8L), N2.plan = c(8L, 8L),
               restricted = FALSE, N.min = NULL, N.max = 60, search = 'crossing',
               search.limit = c(2, 50))
  re <- do.call(strat_reestimate, args)
  hi <- strat_bounds(c(0.5, 0.5), c(1, 1), c(4, 4), c(4, 4), NULL, 60)[2]
  ok <- !is.na(re$n.re)
  expect_true(all(re$n.re[ok] <= hi))
  # Where the normal approximation lies within the bound, the result is that of the
  # search without a bound, capped at the bound
  init <- strat_raw_n(re$hat.p1, re$hat.p2, c(0.5, 0.5), c(1, 1), 0.05, 0.8, 'greater',
                      FALSE, FALSE)
  inside <- ok & ceil_tol(init) <= hi
  expect_gt(sum(inside), 0)
  expect_gt(sum(ok & !inside), 0)
  free <- strat_ss_search(re$hat.p1[inside, , drop = FALSE],
                          re$hat.p2[inside, , drop = FALSE], c(0.5, 0.5), c(1, 1), 0.05,
                          0.8, 'CMH', 'greater', 'minlike')$n
  expect_equal(re$n.re[inside], pmin(free, hi))
})

