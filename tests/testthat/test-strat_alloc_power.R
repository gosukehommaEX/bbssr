p1 <- matrix(c(0.5, 0.8), 1)
p2 <- matrix(c(0.2, 0.4), 1)

test_that("strat_alloc_power reports the probability of the configurations left out", {
  al <- check_strat_alloc('both', 16, c(0.5, 0.5), 1, 4)
  full <- strat_alloc_power(16, al, 'both', p1, p2, 0.05, 'CMH', 'greater', 'minlike',
                            eps = 0)
  cut <- strat_alloc_power(16, al, 'both', p1, p2, 0.05, 'CMH', 'greater', 'minlike',
                           eps = 1e-6)
  expect_identical(full$skipped, 0)
  expect_gt(cut$skipped, 0)
  expect_lt(cut$configs, full$configs)
  expect_lte(abs(full$power - cut$power), cut$skipped)
  expect_null(full$n)
  # 17 stratum sizes, each with (n1 + 1) (n2 + 1) group sizes
  expect_identical(full$configs, as.integer(sum((0:16 + 1) * (16:0 + 1))))
})

test_that("strat_alloc_power keeps the stratum sizes of alloc = 'group'", {
  al <- check_strat_alloc('group', 17, c(0.4, 0.6), 3, 4)
  res <- strat_alloc_power(17, al, 'group', p1, p2, 0.05, 'CMH', 'greater', 'minlike')
  expect_identical(res$n, c(6L, 11L))
  expect_identical(res$configs, 7L * 12L)
  # The rounding tolerates a product slightly below an integer, 100 * 0.29 = 28.999...
  al <- check_strat_alloc('group', 100, c(0.29, 0.71), 1, 4)
  expect_identical(strat_alloc_power(100, al, 'group', p1, p2, 0.05, 'CMH', 'greater',
                                     'minlike', eps = 1e-6)$n, c(29L, 71L))
})

test_that("strat_alloc_power with permuted blocks has few configurations", {
  # With blocks of 4 holding 2 patients of group 1, a stratum of n patients has
  # min(n %% 4, 4 - n %% 4) + 1 possible sizes of group 1
  al <- check_strat_alloc('stratum', 10, c(0.5, 0.5), 1, 4)
  res <- strat_alloc_power(10, al, 'stratum', p1, p2, 0.05, 'CMH', 'greater', 'minlike',
                           eps = 0)
  k <- function(n) pmin(n %% 4, 4 - n %% 4) + 1
  expect_identical(res$configs, as.integer(sum(k(0:10) * k(10:0))))
})
