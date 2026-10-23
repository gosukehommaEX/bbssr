#' Interim Group Sizes of the Strata
#'
#' Internal helper returning the group sizes of the strata at the interim analysis of a
#' stratified design with blinded sample size re-estimation. With \code{omega}, the size
#' of group 2 in stratum \code{j} is \code{ceiling(omega N2[j])} and that of group 1 is
#' \code{ceiling(r[j] n2[j])}, as in the designs without strata. Otherwise the sizes are
#' taken from \code{n.interim}.
#'
#' @param N2 Planned sizes of group 2 in the strata
#' @param r Allocation ratios of the strata
#' @param omega Fraction of the planned sizes at the interim analysis, or \code{NULL}
#' @param n.interim Matrix with two rows, the interim sizes of group 1 and of group 2, and
#'   one column per stratum, or \code{NULL}
#'
#' @return A list with the integer vectors \code{n1} and \code{n2}
#'
#' @keywords internal
#' @noRd
strat_n1 <- function(N2, r, omega, n.interim) {
  J <- length(N2)
  if (!is.null(n.interim)) {
    if (!is.null(omega)) stop('supply either omega or n.interim, not both')
    if (!is.matrix(n.interim) || !is.numeric(n.interim) || nrow(n.interim) != 2 ||
        ncol(n.interim) != J || anyNA(n.interim) || any(n.interim < 1) ||
        any(n.interim != round(n.interim))) {
      stop('n.interim must be a matrix of positive integers with two rows, the interim ',
           'sizes of group 1 and group 2, and one column per stratum')
    }
    return(list(n1 = as.integer(n.interim[1, ]), n2 = as.integer(n.interim[2, ])))
  }
  if (is.null(omega)) stop('supply either omega or n.interim')
  if (length(omega) != 1 || !is.numeric(omega) || is.na(omega) || omega <= 0 ||
      omega > 1) {
    stop('omega must be a single value in (0, 1]')
  }
  n2 <- ceiling(omega * N2)
  list(n1 = as.integer(ceiling(r * n2)), n2 = as.integer(n2))
}
