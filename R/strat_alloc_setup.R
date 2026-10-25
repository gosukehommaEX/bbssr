#' Allocation Sequences of a Re-Estimation Design with Random Sizes
#'
#' Internal helper shared by \code{\link{StratPowerBSSR}}, \code{\link{StratTypeIErrorBSSR}}
#' and \code{\link{StratAlphaAdjBSSR}} with a random allocation and two strata. It checks
#' the settings, draws \code{nseq} allocation sequences of the interim patients, who are
#' the first \code{sum(n1 + n2)} patients of \code{strat_design}, and returns a function
#' that lays out the sequences for a level of the re-estimation.
#'
#' @param d List returned by \code{strat_design}
#' @param alloc \code{'stratum'}, \code{'group'} or \code{'both'}
#' @param block Size of the permuted blocks
#' @param nseq Number of sequences
#' @param re.args Function of the level of the re-estimation returning the list of the
#'   arguments of \code{strat_reestimate}
#' @param ss.frac \code{'observed'} or \code{'planned'}
#' @param restricted Logical. Whether the final total is at least the planned total
#' @param N.min Smallest final total, or \code{NULL}
#' @param N.max Largest final total, or \code{NULL}
#'
#' @return A list with \code{al} of \code{check_strat_alloc}, the planned total
#'   \code{N.plan}, the interim total \code{n} and the function \code{layout} of the level
#'   of the re-estimation, which returns the result of \code{strat_alloc_layout}
#'
#' @keywords internal
#' @noRd
strat_alloc_setup <- function(d, alloc, block, nseq, re.args, ss.frac, restricted, N.min,
                              N.max) {
  if (d$J != 2) stop('a random allocation is available for designs with two strata')
  if (length(nseq) != 1 || !is.numeric(nseq) || is.na(nseq) || nseq < 2 ||
      nseq != round(nseq)) {
    stop('nseq must be a single integer of at least 2')
  }
  al <- check_strat_alloc(alloc, NULL, d$frac, d$r, block)
  a <- d$a
  N.plan <- sum(a$N1 + a$N2)
  n <- sum(d$n1 + d$n2)
  if (!is.null(N.max) && N.max < n) {
    stop('N.max must be at least the interim total sample size')
  }
  seqs <- strat_alloc_seqs(alloc, al, n, nseq)
  layout <- function(level.ss) {
    args <- re.args(level.ss)
    totals <- function(nS) {
      strat_alloc_reestimate(nS, args, ss.frac, N.plan, restricted, N.min, N.max,
                             a$alternative, a$tsmethod)$N
    }
    lay <- strat_alloc_layout(seqs, n, alloc, al, totals)
    # The extended sequences are kept, so that a later level continues the same draws
    seqs <<- lay$seqs
    lay
  }
  list(al = al, N.plan = N.plan, n = n, layout = layout)
}
