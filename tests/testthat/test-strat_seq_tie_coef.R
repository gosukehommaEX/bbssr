# The sequences of test-strat_seq_power.R; reference values from
# tools/reference/reference_strat.py (seq_power with p1 = p2 = theta)
seq_input <- function(sq) {
  code <- 2 * (sq[, 1] - 1) + sq[, 2]
  counts <- rbind(0L, apply(outer(code, 1:4, '=='), 2, cumsum))
  storage.mode(counts) <- 'integer'
  nS <- c(sum(counts[13, 1:2]), sum(counts[13, 3:4]))
  s <- expand.grid(s1 = 0:nS[1], s2 = 0:nS[2])
  list(n = 12L, counts = counts,
       Nfin = as.integer(pmin(14 + 2 * (s$s1 %% 3) + s$s2, 26)))
}
seq_a <- rbind(c(1, 1), c(2, 1), c(1, 2), c(2, 1), c(1, 1), c(1, 2), c(2, 1), c(1, 2),
               c(2, 1), c(1, 1), c(1, 2), c(2, 1), c(2, 2), c(1, 1), c(2, 2), c(1, 2),
               c(2, 2), c(1, 1), c(2, 1), c(2, 2), c(1, 2), c(1, 1), c(2, 2), c(2, 1),
               c(1, 1), c(2, 2))
seq_b <- rbind(c(2, 2), c(1, 1), c(1, 2), c(2, 1), c(1, 1), c(2, 1), c(1, 2), c(1, 2),
               c(2, 2), c(1, 1), c(2, 1), c(1, 2), c(2, 1), c(1, 1), c(2, 2), c(2, 2),
               c(1, 2), c(2, 1), c(1, 1), c(2, 2), c(1, 1), c(2, 1), c(1, 2), c(2, 2),
               c(2, 1), c(1, 2))
seqs <- list(seq_input(seq_a), seq_input(seq_b))

test_that("strat_seq_tie_coef agrees with the reference values given the sequences", {
  for (Test in c('CMH', 'Strat-Fisher')) {
    a <- check_strat_args(c(1, 1), c(1, 1), 0.05, Test, 'greater', 'minlike')
    res <- strat_seq_tie_coef(seqs, a$type, a$alt, a$ts, 0.05, c(0.3, 0.7), c(0.6, 0.2))
    ref <- if (Test == 'CMH') {
      c(0.0519055528697931, 0.0423666910114529, 0.0508571256887412, 0.0455905560693198)
    } else {
      c(0.018306931098221, 0.0119451674496762, 0.018604465633067, 0.011912247243286)
    }
    expect_equal(as.vector(t(res$values)), ref, tolerance = 1e-10)
    # The coefficients average the polynomials of the sequences
    expect_equal(strat_bern_eval(res$coef, c(0.3, 0.7), c(0.6, 0.2)),
                 colMeans(res$values), tolerance = 1e-12)
  }
})

test_that("strat_seq_tie_coef raises the degrees to the largest final stratum totals", {
  res <- strat_seq_tie_coef(seqs, 0L, 0L, 0L, 0.05, numeric(0), numeric(0))
  D <- function(j) {
    max(vapply(seqs, function(s) max(rowSums(s$counts[s$Nfin + 1, 2 * j - 1:0])),
               numeric(1)))
  }
  expect_identical(dim(res$coef), as.integer(c(D(1), D(2)) + 1))
  expect_identical(dim(res$values), c(2L, 0L))
  # A single sequence repeated gives its own polynomial
  one <- strat_seq_tie_coef(seqs[c(1, 1)], 0L, 0L, 0L, 0.05, 0.4, 0.5)
  expect_equal(one$values[1, 1], one$values[2, 1])
})
