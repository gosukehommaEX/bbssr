# Two allocation sequences of 26 patients, as (stratum, group) in the order of enrolment;
# the first has no patient of group 2 in stratum 2 among the first 12. The interim
# analysis has the first 12 patients and the final analysis the first
# min(26, 14 + 2 (s1 %% 3) + s2) patients. Reference values from
# tools/reference/reference_strat.py (seq_power), which applies bssr_power() to the
# group sizes read from the sequences
seq_a <- rbind(c(1, 1), c(2, 1), c(1, 2), c(2, 1), c(1, 1), c(1, 2), c(2, 1), c(1, 2),
               c(2, 1), c(1, 1), c(1, 2), c(2, 1), c(2, 2), c(1, 1), c(2, 2), c(1, 2),
               c(2, 2), c(1, 1), c(2, 1), c(2, 2), c(1, 2), c(1, 1), c(2, 2), c(2, 1),
               c(1, 1), c(2, 2))
seq_b <- rbind(c(2, 2), c(1, 1), c(1, 2), c(2, 1), c(1, 1), c(2, 1), c(1, 2), c(1, 2),
               c(2, 2), c(1, 1), c(2, 1), c(1, 2), c(2, 1), c(1, 1), c(2, 2), c(2, 2),
               c(1, 2), c(2, 1), c(1, 1), c(2, 2), c(1, 1), c(2, 1), c(1, 2), c(2, 2),
               c(2, 1), c(1, 2))
seq_input <- function(sq) {
  code <- 2 * (sq[, 1] - 1) + sq[, 2]
  counts <- rbind(0L, apply(outer(code, 1:4, '=='), 2, cumsum))
  storage.mode(counts) <- 'integer'
  nS <- c(sum(counts[13, 1:2]), sum(counts[13, 3:4]))
  s <- expand.grid(s1 = 0:nS[1], s2 = 0:nS[2])
  list(n = 12L, counts = counts,
       Nfin = as.integer(pmin(14 + 2 * (s$s1 %% 3) + s$s2, 26)))
}
seqs <- list(seq_input(seq_a), seq_input(seq_b))

test_that("strat_seq_power agrees with the reference values given the sequences", {
  pw <- function(Test) {
    a <- check_strat_args(c(1, 1), c(1, 1), 0.05, Test, 'greater', 'minlike')
    res <- strat_seq_power(seqs, rbind(c(0.5, 0.8), c(0.3, 0.6)),
                           rbind(c(0.2, 0.4), c(0.3, 0.6)), a$type, a$alt, a$ts, 0.05)
    expect_identical(dim(res), c(2L, 2L))
    as.vector(t(res))
  }
  expect_equal(pw('CMH'), c(0.496914979099771, 0.0519055528697931, 0.483993049138463,
                            0.0508571256887412), tolerance = 1e-10)
  expect_equal(pw('Strat-Fisher'), c(0.323853159204127, 0.018306931098221,
                                     0.313380191202181, 0.018604465633067),
               tolerance = 1e-10)
  expect_equal(pw('Unstrat-Chisq'), c(0.545591222316565, 0.0674361891351982,
                                      0.492680398724138, 0.050987708156832),
               tolerance = 1e-10)
})

test_that("strat_seq_power checks the sequences", {
  bad <- seqs[[1]]
  bad$counts[5, 1] <- bad$counts[5, 1] + 1L
  expect_error(strat_seq_power(list(bad), matrix(0.5, 1, 2), matrix(0.5, 1, 2), 0, 0, 0,
                               0.05), 'sum to T')
  bad <- seqs[[1]]
  bad$Nfin[1] <- 30L
  expect_error(strat_seq_power(list(bad), matrix(0.5, 1, 2), matrix(0.5, 1, 2), 0, 0, 0,
                               0.05), 'final totals')
  bad <- seqs[[1]]
  bad$Nfin <- bad$Nfin[-1]
  expect_error(strat_seq_power(list(bad), matrix(0.5, 1, 2), matrix(0.5, 1, 2), 0, 0, 0,
                               0.05), 'one element per interim cell')
  expect_error(strat_seq_power(list(), matrix(0.5, 1, 2), matrix(0.5, 1, 2), 0, 0, 0,
                               0.05), 'at least one sequence')
})
