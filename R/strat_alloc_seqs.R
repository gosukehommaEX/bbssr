#' Allocation Sequences of a Two-Stratum Re-Estimation Design with Random Sizes
#'
#' Internal helper drawing \code{nseq} allocation sequences of the first \code{n} patients
#' with \code{strat_alloc_seq}, the patients of the interim analysis. A sequence without
#' interim patients in a stratum is drawn again, so the sequences are those with interim
#' patients in both strata. \code{strat_alloc_layout} extends the sequences to the final
#' totals.
#'
#' @param alloc \code{'stratum'}, \code{'group'} or \code{'both'}
#' @param al List returned by \code{check_strat_alloc}
#' @param n Interim total sample size
#' @param nseq Number of sequences
#'
#' @return A list of \code{nseq} results of \code{strat_alloc_seq}
#'
#' @keywords internal
#' @noRd
strat_alloc_seqs <- function(alloc, al, n, nseq) {
  if (alloc == 'group') {
    k <- floor(n * al$frac[1] + 1e-9)
    if (k == 0 || k == n) {
      stop('the interim analysis of the first ', n, ' patients has no patient in one ',
           'of the strata')
    }
  }
  lapply(seq_len(nseq), function(i) {
    for (try in 1:1000) {
      x <- strat_alloc_seq(alloc, al, n)
      c0 <- x$counts[n + 1, ]
      if (c0[1] + c0[2] > 0 && c0[3] + c0[4] > 0) return(x)
    }
    stop('the interim analyses of 1000 sequences had no patient in one of the strata')
  })
}
