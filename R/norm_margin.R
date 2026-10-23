#' Margins of a Test of Equivalence
#'
#' Internal helper returning the margin in the form used by the computation. For
#' \code{alternative = 'equivalence'} it validates the margin and returns the pair
#' \code{c(lower, upper)}. On the scale of the risk difference the null hypothesis is
#' \code{p1 - p2 <= -lower} or \code{p1 - p2 >= upper}, both margins lie in (0, 1), and a
#' single value \code{m} gives \code{c(m, m)}. On the scale of the risk ratio the null
#' hypothesis is \code{p1 / p2 <= lower} or \code{p1 / p2 >= upper} with
#' \code{0 < lower < 1 < upper}, and a single ratio \code{R} gives
#' \code{c(min(R, 1 / R), max(R, 1 / R))}. For the other alternatives the margin is
#' returned unchanged; it is validated by \code{check_rr_args} and \code{check_delta}.
#'
#' @param alternative \code{'greater'}, \code{'less'}, \code{'two.sided'} or
#'   \code{'equivalence'}
#' @param margin Margin as supplied to a function of the package
#' @param margin.scale \code{'RD'} or \code{'RR'}
#'
#' @return The margin, a numeric vector of length 2 for a test of equivalence
#'
#' @keywords internal
#' @noRd
norm_margin <- function(alternative, margin, margin.scale) {
  if (!identical(alternative, 'equivalence')) return(margin)
  if (!is.numeric(margin) || !(length(margin) %in% 1:2) || anyNA(margin) ||
      any(!is.finite(margin))) {
    stop('the margin of a test of equivalence must be one or two finite numbers')
  }
  if (identical(margin.scale, 'RR')) {
    if (length(margin) == 1) {
      if (margin <= 0 || margin == 1) {
        stop('a single margin on the scale of the risk ratio must be positive and ',
             'differ from 1')
      }
      margin <- c(min(margin, 1 / margin), max(margin, 1 / margin))
    }
    if (!(margin[1] > 0 && margin[1] < 1 && margin[2] > 1)) {
      stop('the margins c(lower, upper) on the scale of the risk ratio must satisfy ',
           '0 < lower < 1 < upper')
    }
  } else {
    if (length(margin) == 1) margin <- c(margin, margin)
    if (any(margin <= 0) || any(margin >= 1)) {
      stop('the margins c(lower, upper) on the scale of the risk difference must lie ',
           'in (0, 1)')
    }
  }
  as.double(margin)
}
