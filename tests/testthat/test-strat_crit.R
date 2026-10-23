test_that("strat_crit agrees with a scan of the p-values of every outcome", {
  N1 <- c(4L, 6L)
  N2 <- c(5L, 3L)
  z <- as.matrix(expand.grid(0:9, 0:9))
  # One table for every pair (z, S): the responders of group 1 fill stratum 1 first
  rows <- do.call(rbind, lapply(seq_len(nrow(z)), function(k) {
    lo <- pmax(0, z[k, ] - N2)
    hi <- pmin(z[k, ], N1)
    t(vapply(sum(lo):sum(hi), function(S) {
      add1 <- min(S - sum(lo), hi[1] - lo[1])
      c(k, S, lo + c(add1, S - sum(lo) - add1))
    }, numeric(4)))
  }))
  n <- nrow(rows)
  x1 <- matrix(as.integer(rows[, 3:4]), n)
  x2 <- matrix(as.integer(z[rows[, 1], ] - rows[, 3:4]), n)
  for (type in 0:4) {
    for (alt in 0:2) {
      ts <- if (type == 4 && alt == 2) 2L else 0L
      crit <- strat_crit(N1, N2, type, alt, ts, 0.1)
      expect_length(crit$lo, nrow(z))
      p <- strat_pvalue_cases(matrix(N1, n, 2, byrow = TRUE), matrix(N2, n, 2, byrow = TRUE),
                              x1, x2, type, alt, ts)
      inside <- rows[, 2] <= crit$lo[rows[, 1]] | rows[, 2] >= crit$hi[rows[, 1]]
      expect_identical(p %<<% 0.1, inside, info = paste(type, alt))
      expect_gt(sum(inside), 0)
    }
  }
})

test_that("strat_crit checks its arguments", {
  expect_error(strat_crit(c(4L, 6L), c(5L, 0L), 0L, 0L, 0L, 0.1), 'positive')
  expect_error(strat_crit(c(4L, 6L), c(5L, 3L), 5L, 0L, 0L, 0.1), 'type')
  expect_error(strat_crit(c(4L, 6L), c(5L, 3L), 0L, 3L, 0L, 0.1), 'alternative')
  expect_error(strat_crit(c(4L, 6L), c(5L, 3L), 0L, 0L, 0L, 1), 'alpha')
})
