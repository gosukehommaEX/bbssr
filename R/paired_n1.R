#' Number of Pairs at the Interim Analysis
#'
#' Internal helper returning the number of pairs at the interim analysis of a design with
#' blinded sample size re-estimation for paired binary data, given either the fraction
#' \code{omega} of the planned number of pairs \code{N}, rounded up, or the number
#' \code{n.interim} itself.
#'
#' @param N Planned number of pairs, used with \code{omega}
#' @param omega Fraction of \code{N} at the interim analysis, or \code{NULL}
#' @param n.interim Number of pairs at the interim analysis, or \code{NULL}
#'
#' @return A positive integer
#'
#' @keywords internal
#' @noRd
paired_n1 <- function(N, omega, n.interim) {
  if (!is.null(n.interim)) {
    if (!is.null(omega)) stop('supply either omega or n.interim, not both')
    if (length(n.interim) != 1 || is.na(n.interim) || n.interim < 1 ||
        n.interim != round(n.interim)) {
      stop('n.interim must be a single positive integer, the number of interim pairs')
    }
    return(as.integer(n.interim))
  }
  if (is.null(omega)) stop('supply either omega or n.interim')
  if (length(omega) != 1 || is.na(omega) || omega <= 0 || omega > 1) {
    stop('omega must be a single value in (0, 1]')
  }
  if (is.null(N) || length(N) != 1 || is.na(N) || N < 1 || N != round(N)) {
    stop('N must be a single positive integer when omega is given')
  }
  as.integer(ceiling(omega * N))
}
