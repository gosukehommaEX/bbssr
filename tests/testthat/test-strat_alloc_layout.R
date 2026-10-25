test_that("strat_alloc_layout extends the sequences to the final totals of their cells", {
  set.seed(9)
  al <- check_strat_alloc('stratum', NULL, c(0.5, 0.5), 1, 4)
  seqs <- strat_alloc_seqs('stratum', al, 10, 3)
  calls <- 0
  totals <- function(nS) {
    calls <<- calls + 1
    as.integer(10 + seq_len((nS[1] + 1) * (nS[2] + 1)) %% 7)
  }
  lay <- strat_alloc_layout(seqs, 10, 'stratum', al, totals)
  nS <- t(vapply(seqs, function(s) c(sum(s$counts[11, 1:2]), sum(s$counts[11, 3:4])),
                 numeric(2)))
  # The totals are computed once for each pair of interim stratum sizes
  expect_equal(calls, nrow(unique(nS)))
  for (i in 1:3) {
    s <- lay$inputs[[i]]
    expect_identical(s$n, 10L)
    expect_identical(s$Nfin, totals(nS[i, ]))
    expect_identical(nrow(s$counts), max(s$Nfin) + 1L)
    expect_identical(s$counts[1:11, ], seqs[[i]]$counts)
    expect_identical(s$n1, s$counts[11, c(1, 3)])
    expect_identical(s$n2, s$counts[11, c(2, 4)])
    expect_identical(lay$seqs[[i]]$counts, s$counts)
  }
})
