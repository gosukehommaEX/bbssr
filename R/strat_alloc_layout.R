#' Final Totals of the Allocation Sequences of a Re-Estimation Design
#'
#' Internal helper giving, for each sequence of \code{strat_alloc_seqs}, the final totals
#' of the interim cells from the interim stratum sizes, computed once for each pair of
#' sizes, and extending the sequence to the largest of them.
#'
#' @param seqs List of results of \code{strat_alloc_seq}, the interim patients first
#' @param n Interim total sample size
#' @param alloc \code{'stratum'}, \code{'group'} or \code{'both'}
#' @param al List returned by \code{check_strat_alloc}
#' @param totals Function of the interim stratum sizes returning the final totals of the
#'   interim cells, \code{s1} varying fastest
#'
#' @return A list with the extended sequences \code{seqs} and the list \code{inputs} of
#'   the sequences as \code{strat_seq_power} expects them: \code{n}, the integer matrix
#'   \code{counts} and the integer vector \code{Nfin}, with the interim group sizes
#'   \code{n1} and \code{n2}
#'
#' @keywords internal
#' @noRd
strat_alloc_layout <- function(seqs, n, alloc, al, totals) {
  cache <- new.env()
  inputs <- vector('list', length(seqs))
  for (i in seq_along(seqs)) {
    c0 <- seqs[[i]]$counts[n + 1, ]
    nS <- c(c0[1] + c0[2], c0[3] + c0[4])
    key <- paste(nS, collapse = '_')
    if (is.null(cache[[key]])) cache[[key]] <- totals(nS)
    Nfin <- cache[[key]]
    seqs[[i]] <- strat_alloc_seq(alloc, al, max(Nfin), seqs[[i]])
    rows <- seq_len(max(Nfin) + 1)
    inputs[[i]] <- list(n = as.integer(n), counts = seqs[[i]]$counts[rows, , drop = FALSE],
                        Nfin = Nfin, n1 = c0[c(1, 3)], n2 = c0[c(2, 4)])
  }
  list(seqs = seqs, inputs = inputs)
}
