test_that("strat_ss_search finds the crossing of the target power", {
  p1 <- rbind(c(0.65, 0.8))
  p2 <- rbind(c(0.3, 0.45))
  s <- strat_ss_search(p1, p2, c(0.5, 0.5), c(1, 1), 0.025, 0.8, 'CMH', 'greater',
                       'minlike')
  pw <- function(n) {
    sz <- strat_sizes(n, c(0.5, 0.5), c(1, 1))
    strat_fixed_power(sz$N1[1, ], sz$N2[1, ], p1, p2, 0.025, 'CMH', 'greater', 'minlike')
  }
  expect_gte(pw(s$n), 0.8)
  expect_lt(pw(s$n - 1), 0.8)
  expect_true(is.na(s$limit))
})

test_that("strat_ss_search supports the stable search", {
  p1 <- rbind(c(0.65, 0.8))
  p2 <- rbind(c(0.3, 0.45))
  st <- strat_ss_search(p1, p2, c(0.5, 0.5), c(1, 1), 0.025, 0.8, 'CMH', 'greater',
                        'minlike', search = 'stable', search.limit = c(1.2, 10))
  pw <- function(n) {
    sz <- strat_sizes(n, c(0.5, 0.5), c(1, 1))
    strat_fixed_power(sz$N1[1, ], sz$N2[1, ], p1, p2, 0.025, 'CMH', 'greater', 'minlike')
  }
  expect_false(is.na(st$limit))
  expect_true(all(vapply(st$n:st$limit, pw, numeric(1)) >= 0.8 - 1e-8))
  expect_lt(pw(st$n - 1), 0.8)
  expect_error(strat_ss_search(p1, p2, c(0.5, 0.5), c(1, 1), 0.025, 0.8, 'CMH', 'greater',
                               'minlike', search = 'stable', search.limit = c(0.5, 10)),
               'search.limit')
  expect_error(strat_ss_search(p2, p2, c(0.5, 0.5), c(1, 1), 0.025, 0.8, 'CMH', 'greater',
                               'minlike'), 'no effect')
})
