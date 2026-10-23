test_that("paired_ss_search finds the crossing of the target power", {
  # Two pairs of probabilities share the rejection regions of one search
  s <- paired_ss_search(c(0.3, 0.35), c(0.1, 0.1), 0.025, 0.8, 'McNemar-exact',
                        'greater', 100, 0, FALSE, 0)
  expect_equal(s$N, c(84, 60))
  expect_true(all(is.na(s$limit)))
  pw <- function(N, p10, p01, ...) {
    paired_power(get_paired_rr(N, 0.025, ..., 'greater', 100, 0, FALSE, 0), p10, p01)
  }
  expect_gte(pw(84, 0.3, 0.1, 'McNemar-exact'), 0.8)
  expect_lt(pw(83, 0.3, 0.1, 'McNemar-exact'), 0.8)
  expect_gte(pw(60, 0.35, 0.1, 'McNemar-exact'), 0.8)
  expect_lt(pw(59, 0.35, 0.1, 'McNemar-exact'), 0.8)
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
