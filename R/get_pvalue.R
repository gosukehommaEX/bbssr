#' p-Values Shared Across Calls
#'
#' Internal helper returning the matrix of p-values computed by \code{rr_pvalue}. The
#' p-values depend only on the sample sizes and the test, and not on the level, the
#' response probabilities, the interim fraction or the assumed treatment effect. Sample
#' size searches, re-estimation designs and searches for an adjusted level request the
#' same matrices many times, so each matrix is computed once and kept for the rest of the
#' session. A rejection region at any level then follows by comparison.
#'
#' The matrices are held in the environment \code{.bbssr_cache}, keyed by every argument.
#' The store is emptied when it would exceed \code{.bbssr_cache$max.cells} cells, and
#' \code{options(bbssr.cache = FALSE)} disables it.
#'
#' @param N1 Sample size for group 1, an integer
#' @param N2 Sample size for group 2, an integer
#' @param Test Full name of the test
#' @param alternative \code{'greater'}, \code{'less'} or \code{'two.sided'}
#' @param tsmethod \code{'minlike'} or \code{'central'}
#' @param n.grid Number of grid points over the nuisance parameter, an integer
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure
#'
#' @return A numeric matrix of dimension \code{(N1 + 1)} by \code{(N2 + 1)}
#'
#' @keywords internal
#' @noRd
get_pvalue <- function(N1, N2, Test, alternative, tsmethod, n.grid, bb.gamma) {
  if (!isTRUE(getOption('bbssr.cache', TRUE))) {
    return(rr_pvalue(N1, N2, Test, alternative, tsmethod, n.grid, bb.gamma))
  }
  num <- function(x) sprintf('%.17g', as.double(x))
  key <- paste(num(N1), num(N2), Test, alternative, tsmethod, num(n.grid), num(bb.gamma),
               sep = '|')
  hit <- .bbssr_cache$pv[[key]]
  if (!is.null(hit)) return(hit)
  pv <- rr_pvalue(N1, N2, Test, alternative, tsmethod, n.grid, bb.gamma)
  if (.bbssr_cache$cells + length(pv) > .bbssr_cache$max.cells) clear_rr_cache()
  assign(key, pv, envir = .bbssr_cache$pv)
  .bbssr_cache$cells <- .bbssr_cache$cells + length(pv)
  pv
}
