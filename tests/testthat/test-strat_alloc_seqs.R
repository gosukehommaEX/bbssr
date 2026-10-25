test_that("strat_alloc_seqs draws sequences with interim patients in both strata", {
  # With 5 interim patients and a stratum of probability 0.1, most first draws leave
  # stratum 1 empty and are drawn again
  set.seed(8)
  al <- check_strat_alloc('both', NULL, c(0.1, 0.9), 1, 4)
  x <- strat_alloc_seqs('both', al, 5, 50)
  expect_length(x, 50)
  ok <- vapply(x, function(s) {
    c0 <- s$counts[6, ]
    nrow(s$counts) == 6 && c0[1] + c0[2] > 0 && c0[3] + c0[4] > 0
  }, logical(1))
  expect_true(all(ok))
})

test_that("strat_alloc_seqs stops when the fixed stratum sizes leave a stratum empty", {
  al <- check_strat_alloc('group', NULL, c(0.05, 0.95), 1, 4)
  expect_error(strat_alloc_seqs('group', al, 10, 2), 'no patient in one of the strata')
  x <- strat_alloc_seqs('group', al, 20, 2)
  expect_identical(x[[1]]$counts[21, 1] + x[[1]]$counts[21, 2], 1L)
})
