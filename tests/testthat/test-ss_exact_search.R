test_that("ss_exact_search reproduces separate searches for every pair", {
  p1 <- c(0.5, 0.6, 0.45, 0.7)
  p2 <- c(0.2, 0.3, 0.25, 0.3)
  for (tst in c('Chisq', 'Fisher', 'Boschloo')) {
    res <- ss_exact_search(p1, p2, 1, 0.025, 0.8, tst, 'greater', 'minlike', 100L, 0,
                           FALSE, 0)
    got <- res$N2
    expect_type(got, 'integer')
    expect_true(all(is.na(res$limit)))
    one <- vapply(seq_along(p1), function(k) {
      ss_exact_search(p1[k], p2[k], 1, 0.025, 0.8, tst, 'greater', 'minlike', 100L, 0,
                      FALSE, 0)$N2
    }, integer(1))
    expect_identical(got, one, info = tst)
  }
})

test_that("the size attains the target power and one patient fewer does not", {
  N2 <- ss_exact_search(0.6, 0.25, 2, 0.025, 0.8, 'Chisq', 'greater', 'minlike', 100L, 0,
                        FALSE, 0)$N2
  at <- BinaryPower(0.6, 0.25, ceiling(2 * N2), N2, 0.025, 'Chisq')$Power
  below <- BinaryPower(0.6, 0.25, ceiling(2 * (N2 - 1)), N2 - 1, 0.025, 'Chisq')$Power
  expect_gte(at, 0.8)
  expect_lt(below, 0.8)
  expect_identical(ss_exact_search(numeric(0), numeric(0), 1, 0.025, 0.8, 'Chisq',
                                   'greater', 'minlike', 100L, 0, FALSE, 0)$N2,
                   integer(0))
})

test_that("the three searches follow their definitions", {
  run <- function(p1, p2, r, tar.power, Test, search) {
    ss_exact_search(p1, p2, r, 0.025, tar.power, Test, 'greater', 'minlike', 100L, 0,
                    FALSE, 0, search)
  }
  # Reference values from tools/reference/reference_values.py: the sizes of group 2 found
  # by the three searches and the limit of the stable search
  cr <- run(0.6, 0.2, 2, 0.9, 'Chisq', 'crossing')
  sm <- run(0.6, 0.2, 2, 0.9, 'Chisq', 'smallest')
  st <- run(0.6, 0.2, 2, 0.9, 'Chisq', 'stable')
  expect_equal(c(cr$N2, sm$N2, st$N2, st$limit), c(23, 21, 23, 72))
  expect_true(is.na(cr$limit) && is.na(sm$limit))
  # The smallest size attains the target power and no smaller size does
  pw <- vapply(1:21, function(n) {
    BinaryPower(0.6, 0.2, 2 * n, n, 0.025, 'Chisq')$Power
  }, numeric(1))
  expect_true(all(pw[1:20] < 0.9))
  expect_gte(pw[21], 0.9)
  cr <- run(0.6, 0.3, 1, 0.85, 'Fisher', 'crossing')
  sm <- run(0.6, 0.3, 1, 0.85, 'Fisher', 'smallest')
  st <- run(0.6, 0.3, 1, 0.85, 'Fisher', 'stable')
  expect_equal(c(cr$N2, sm$N2, st$N2, st$limit), c(52, 52, 56, 98))
  # Every size from the stable size to the limit attains the target power and the size
  # one unit smaller does not
  pw <- vapply(55:98, function(n) {
    BinaryPower(0.6, 0.3, n, n, 0.025, 'Fisher')$Power
  }, numeric(1))
  expect_lt(pw[1], 0.85)
  expect_true(all(pw[-1] >= 0.85))
})

test_that("the stable search checks its limit", {
  run <- function(limit) {
    ss_exact_search(0.6, 0.2, 2, 0.025, 0.9, 'Chisq', 'greater', 'minlike', 100L, 0,
                    FALSE, 0, 'stable', limit)
  }
  # The size 22 from the normal approximation falls short of the target power
  expect_error(run(c(1, 0)), 'raise search.limit')
  for (bad in list(c(0.5, 10), c(2, -1), 2, c(2, NA), c('2', '50'))) {
    expect_error(run(bad), 'search.limit must be')
  }
})
