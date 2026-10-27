#' Critical Values of a Stratified Test, Stored for Reuse
#'
#' Internal helper returning the critical values of the test of a stratified design for
#' every vector of stratum totals, computed by \code{strat_crit} and stored for the rest
#' of the session like the p-values of the other designs (see \code{\link{bbssr}}). Given
#' the numbers \code{z_j} of responders of the strata, the test rejects when
#' \code{S = sum_j x_j}, the number of responders of group 1, is at most \code{lo} or at
#' least \code{hi}.
#'
#' @inheritParams StratRR
#'
#' @return A list with the integer vectors \code{N1} and \code{N2} and the integer vectors
#'   \code{lo} and \code{hi}, with one element for each vector of stratum totals, the
#'   first stratum varying fastest
#'
#' @keywords internal
#' @noRd
get_strat_crit <- function(N1, N2, alpha, Test, alternative, tsmethod) {
  a <- check_strat_args(N1, N2, alpha, Test, alternative, tsmethod)
  cells <- prod(as.double(a$N1 + a$N2 + 1))
  if (cells > 2e7) {
    stop('the design has too many vectors of stratum totals (', format(cells),
         ') for its critical values to be stored; in a re-estimation design, bound the ',
         "final total with N.max or use ss.method = 'standard'")
  }
  compute <- function() {
    crit <- strat_crit(a$N1, a$N2, a$type, a$alt, a$ts, alpha)
    list(N1 = a$N1, N2 = a$N2, lo = crit$lo, hi = crit$hi)
  }
  if (!isTRUE(getOption('bbssr.cache', TRUE))) return(compute())
  key <- paste('strat', paste(a$N1, collapse = ','), paste(a$N2, collapse = ','), a$Test,
               a$alternative, a$tsmethod, sprintf('%.17g', alpha), sep = '|')
  hit <- .bbssr_cache$pv[[key]]
  if (!is.null(hit)) return(hit)
  out <- compute()
  if (.bbssr_cache$cells + 2 * cells > .bbssr_cache$max.cells) clear_rr_cache()
  assign(key, out, envir = .bbssr_cache$pv)
  .bbssr_cache$cells <- .bbssr_cache$cells + 2 * cells
  out
}
