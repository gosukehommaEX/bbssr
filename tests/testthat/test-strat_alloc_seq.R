test_that("strat_alloc_seq counts the patients of the four groups", {
  set.seed(1)
  al <- check_strat_alloc('stratum', NULL, c(0.4, 0.6), 1, 4)
  x <- strat_alloc_seq('stratum', al, 50)
  cnt <- x$counts
  expect_identical(dim(cnt), c(51L, 4L))
  expect_identical(as.integer(rowSums(cnt)), 0:50)
  expect_true(all(diff(cnt) >= 0))
  expect_identical(cnt[51, ], as.integer(table(factor(2 * (x$st - 1) + x$grp, 1:4))))
})

test_that("strat_alloc_seq allocates permuted blocks within the strata", {
  set.seed(2)
  al <- check_strat_alloc('stratum', NULL, c(0.5, 0.5), 3, 4)
  x <- strat_alloc_seq('stratum', al, 80)
  for (j in 1:2) {
    g <- x$grp[x$st == j]
    k <- floor(length(g) / 4)
    # Every complete block holds 3 patients of group 1 and 1 of group 2
    blocks <- matrix(g[seq_len(4 * k)], nrow = 4)
    expect_true(all(colSums(blocks == 1) == 3))
  }
})

test_that("strat_alloc_seq keeps the stratum sizes of alloc = 'group'", {
  set.seed(3)
  al <- check_strat_alloc('group', NULL, c(0.3, 0.7), 1, 4)
  x <- strat_alloc_seq('group', al, 40)
  T <- 0:40
  expect_identical(x$counts[, 1] + x$counts[, 2], as.integer(floor(T * 0.3 + 1e-9)))
})

test_that("strat_alloc_seq extends a sequence without changing its beginning", {
  al <- check_strat_alloc('both', NULL, c(0.5, 0.5), 1, 4)
  set.seed(4)
  x <- strat_alloc_seq('both', al, 20)
  y <- strat_alloc_seq('both', al, 45, x)
  expect_identical(y$counts[1:21, ], x$counts)
  expect_identical(dim(y$counts), c(46L, 4L))
  set.seed(4)
  expect_identical(strat_alloc_seq('both', al, 20)$counts, x$counts)
})
