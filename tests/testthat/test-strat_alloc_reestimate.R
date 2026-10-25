re.args <- list(Delta.A = 0.3, effect = 'RD', frac = c(0.5, 0.5), r = c(1, 1),
                ss.method = 'standard', ss.Test = 'CMH', ss.alpha = 0.05, tar.power = 0.8,
                N1.plan = c(6, 6), N2.plan = c(6, 6), restricted = FALSE, N.min = NULL,
                N.max = NULL, search = 'crossing', search.limit = c(2, 50),
                ss.strata = 'stratified')

test_that("strat_alloc_reestimate takes the total of the re-estimated group sizes", {
  nS <- c(5, 7)
  res <- strat_alloc_reestimate(nS, re.args, 'observed', 24, FALSE, NULL, NULL, 'greater',
                                'minlike')
  s <- as.matrix(expand.grid(0:5, 0:7))
  expect_length(res$N, 6 * 8)
  re <- strat_reestimate(s, c(3, 4), c(2, 3), 0.3, 'RD', nS / 12, c(1, 1), 'standard',
                         'CMH', 0.05, 0.8, 'greater', 'minlike', c(3, 4), c(2, 3), FALSE,
                         NULL, NULL, 'crossing', c(2, 50))
  ok <- !is.na(re$n.re)
  sz <- strat_sizes(re$n.re[ok], nS / 12, c(1, 1))
  expect_identical(res$N[ok], as.integer(pmax(rowSums(sz$N1) + rowSums(sz$N2), 12)))
  # The planned fractions instead of the observed ones
  pl <- strat_alloc_reestimate(nS, re.args, 'planned', 24, FALSE, NULL, NULL, 'greater',
                               'minlike')
  expect_false(identical(pl$N, res$N))
})

test_that("strat_alloc_reestimate applies the bounds to the total", {
  nS <- c(6, 6)
  free <- strat_alloc_reestimate(nS, re.args, 'observed', 24, FALSE, NULL, NULL, 'greater',
                                 'minlike')$N
  expect_true(all(free >= 12))
  capped <- strat_alloc_reestimate(nS, re.args, 'observed', 24, FALSE, NULL, 30,
                                   'greater', 'minlike')$N
  expect_identical(capped, pmin(free, 30L))
  raised <- strat_alloc_reestimate(nS, re.args, 'observed', 24, TRUE, 28, NULL, 'greater',
                                   'minlike')$N
  expect_identical(raised, pmax(free, 28L))
})
