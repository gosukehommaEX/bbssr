test_that("paired_ss_search finds the crossing of the target power", {
  s <- paired_ss_search(c(0.3, 0.1), c(0.1, 0.1), 0.025, 0.8,
                        c('McNemar-exact'), 'greater', 100, 0, FALSE, 0)
  expect_equal(s$N[1], 84)
  expect_true(all(is.na(s$limit)))
  pw <- function(N, p10, p01, ...) {
    paired_power(get_paired_rr(N, 0.025, ..., 'greater', 100, 0, FALSE, 0), p10, p01)
  }
  expect_gte(pw(84, 0.3, 0.1, 'McNemar-exact'), 0.8)
  expect_lt(pw(83, 0.3, 0.1, 'McNemar-exact'), 0.8)
})

test_that("paired_ss_search supports the three searches", {
  sm <- paired_ss_search(0.3, 0.1, 0.025, 0.8, 'McNemar-midP', 'greater', 100, 0, FALSE, 0,
                         search = 'smallest')
  expect_equal(sm$N, 76)
  st <- paired_ss_search(0.3, 0.1, 0.025, 0.8, 'McNemar-midP', 'greater', 100, 0, FALSE, 0,
                         search = 'stable', search.limit = c(1.5, 10))
  expect_false(is.na(st$limit))
  rr <- function(N) get_paired_rr(N, 0.025, 'McNemar-midP', 'greater', 100, 0, FALSE, 0)
  # Every number from the result up to the limit attains the target power
  pw <- vapply(st$N:st$limit, function(N) paired_power(rr(N), 0.3, 0.1), numeric(1))
  expect_true(all(pw >= 0.8 - 1e-8))
  expect_lt(paired_power(rr(st$N - 1), 0.3, 0.1), 0.8)
  expect_error(paired_ss_search(0.3, 0.1, 0.025, 0.8, 'McNemar', 'greater', 100, 0, FALSE,
                                0, search = 'stable', search.limit = c(0.5, 10)),
               'search.limit')
})
