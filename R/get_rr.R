#' Rejection Region Shared Across Calls
#'
#' Internal helper returning the rejection region of \code{\link{BinaryRR}} as a plain
#' logical matrix. A rejection region depends only on the sample sizes, the level and the
#' test, and not on the response probabilities, the interim fraction or the assumed
#' treatment effect. Sample size searches and re-estimation designs request the same
#' regions many times, so each region is computed once and kept for the rest of the
#' session.
#'
#' The regions are held in the environment \code{.bbssr_cache}, keyed by every argument
#' of \code{\link{BinaryRR}}. The store is emptied when it would exceed
#' \code{.bbssr_cache$max.cells} cells. Setting \code{options(bbssr.cache = FALSE)}
#' disables the store, in which case every call computes the region afresh.
#'
#' @param N1 Sample size for group 1
#' @param N2 Sample size for group 2
#' @param alpha Level of significance
#' @param Test Type of statistical test
#' @param alternative Direction of the alternative hypothesis
#' @param tsmethod Convention used to construct the two-sided conditional tests
#' @param n.grid Number of grid points over the nuisance parameter
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure
#'
#' @return A logical matrix of dimension \code{(N1 + 1)} by \code{(N2 + 1)} without
#'   attributes other than \code{dim}
#'
#' @keywords internal
#' @noRd
get_rr <- function(N1, N2, alpha, Test, alternative, tsmethod, n.grid, bb.gamma) {
  Test <- match.arg(Test, c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo'))
  args <- list(N1, N2, alpha, alternative, tsmethod, n.grid, bb.gamma)
  scalar <- all(vapply(args, function(a) length(a) == 1L && !is.na(a), logical(1)))
  use.cache <- scalar && isTRUE(getOption('bbssr.cache', TRUE))
  if (use.cache) {
    num <- function(x) sprintf('%.17g', as.double(x))
    key <- paste(num(N1), num(N2), num(alpha), Test, alternative, tsmethod,
                 num(n.grid), num(bb.gamma), sep = '|')
    hit <- .bbssr_cache$rr[[key]]
    if (!is.null(hit)) {
      # The warning of BinaryRR() is repeated, so a stored region behaves like a new one
      if (bb.gamma > 0 && !(Test %in% c('Z-pool', 'Boschloo'))) {
        warning('bb.gamma applies only to the Z-pool and Boschloo tests and is ignored here')
      }
      return(hit)
    }
  }
  RR <- BinaryRR(N1, N2, alpha, Test, alternative, tsmethod, n.grid, bb.gamma)
  rr <- matrix(as.vector(RR), nrow = attr(RR, 'N1') + 1L, ncol = attr(RR, 'N2') + 1L)
  if (use.cache) {
    if (.bbssr_cache$cells + length(rr) > .bbssr_cache$max.cells) clear_rr_cache()
    assign(key, rr, envir = .bbssr_cache$rr)
    .bbssr_cache$cells <- .bbssr_cache$cells + length(rr)
  }
  rr
}
