#' Cached p-Values of a Test for Paired Binary Data
#'
#' Internal helper returning \code{paired_pvalue} for the given arguments. The p-values
#' are kept in the cache of the package, which \code{get_pvalue} shares, under a key that
#' starts with \code{'paired'}, so that a later call with the same arguments reuses them.
#' Setting \code{options(bbssr.cache = FALSE)} disables the cache.
#'
#' @inheritParams PairedRR
#'
#' @return The matrix of p-values of \code{paired_pvalue}
#'
#' @keywords internal
#' @noRd
get_paired_pvalue <- function(N, Test, alternative, n.grid, bb.gamma, ref.pvalue, margin) {
  if (!isTRUE(getOption('bbssr.cache', TRUE))) {
    return(paired_pvalue(N, Test, alternative, n.grid, bb.gamma, ref.pvalue, margin))
  }
  num <- function(x) sprintf('%.17g', as.double(x))
  key <- paste('paired', num(N), Test, alternative, num(n.grid), num(bb.gamma), ref.pvalue,
               paste(num(margin), collapse = ','), sep = '|')
  hit <- .bbssr_cache$pv[[key]]
  if (!is.null(hit)) return(hit)
  pv <- paired_pvalue(N, Test, alternative, n.grid, bb.gamma, ref.pvalue, margin)
  if (.bbssr_cache$cells + length(pv) > .bbssr_cache$max.cells) clear_rr_cache()
  assign(key, pv, envir = .bbssr_cache$pv)
  .bbssr_cache$cells <- .bbssr_cache$cells + length(pv)
  pv
}
