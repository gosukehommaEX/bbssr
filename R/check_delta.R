#' Validate an Assumed Treatment Effect
#'
#' Internal helper checking that an assumed treatment effect is a single value on the
#' scale given by \code{effect}, that it differs from the value under the null hypothesis
#' and that its direction agrees with \code{alternative}. With a non-inferiority margin
#' the null value is \code{-margin} for \code{'greater'} and \code{margin} for
#' \code{'less'}, and only the risk difference and one-sided alternatives are allowed.
#' With a margin on the scale of the risk ratio the null value is \code{margin} for both
#' one-sided alternatives, and only the risk ratio is allowed. For a test of equivalence
#' the effect must lie strictly inside the equivalence region, \code{(-lower, upper)} on
#' the scale of the risk difference and \code{(lower, upper)} on that of the risk ratio,
#' and must be on the scale of the margins.
#'
#' @param Delta Treatment effect
#' @param effect \code{'RD'}, \code{'RR'} or \code{'OR'}
#' @param alternative \code{'greater'}, \code{'less'}, \code{'two.sided'} or
#'   \code{'equivalence'}
#' @param name Name of the argument, used in the error messages
#' @param margin Non-inferiority margin, 0 for a test of superiority on the scale of the
#'   risk difference
#' @param margin.scale \code{'RD'} or \code{'RR'}, the scale of the margin
#'
#' @return \code{NULL}, invisibly
#'
#' @keywords internal
#' @noRd
check_delta <- function(Delta, effect, alternative, name = 'Delta.A', margin,
                        margin.scale) {
  if (length(Delta) != 1 || is.na(Delta)) stop(name, ' must be a single value')
  if (alternative == 'equivalence') {
    m <- norm_margin(alternative, margin, margin.scale)
    if (margin.scale == 'RR') {
      if (effect != 'RR') stop("margin.scale = 'RR' requires effect = 'RR'")
      if (!(Delta > m[1] && Delta < m[2])) {
        stop(name, ' must lie strictly between the two margins of the risk ratio')
      }
    } else {
      if (effect != 'RD') {
        stop("a test of equivalence on the scale of the risk difference requires ",
             "effect = 'RD'")
      }
      if (!(Delta > -m[1] && Delta < m[2])) {
        stop(name, ' must lie strictly between -lower and upper, the margins of ',
             'equivalence')
      }
    }
    return(invisible(NULL))
  }
  if (margin.scale == 'RR') {
    if (length(margin) != 1 || !is.numeric(margin) || !is.finite(margin) || margin <= 0) {
      stop("margin must be a single positive value when margin.scale is 'RR'")
    }
    if (effect != 'RR') stop("margin.scale = 'RR' requires effect = 'RR'")
    if (alternative == 'two.sided') {
      stop("margin.scale = 'RR' requires a one-sided alternative")
    }
    if (Delta <= 0) stop(name, ' must be positive for a risk ratio')
    if (alternative == 'greater' && !(Delta > margin)) {
      stop(name, " must exceed margin when alternative is 'greater'")
    }
    if (alternative == 'less' && !(Delta < margin)) {
      stop(name, " must fall below margin when alternative is 'less'")
    }
    return(invisible(NULL))
  }
  if (length(margin) != 1 || !is.numeric(margin) || is.na(margin) || abs(margin) >= 1) {
    stop('margin must be a single value in (-1, 1)')
  }
  if (margin != 0) {
    if (effect != 'RD') stop("a non-zero margin requires effect = 'RD'")
    if (alternative == 'two.sided') {
      stop('a non-zero margin requires a one-sided alternative')
    }
    if (abs(Delta) >= 1) stop(name, ' must lie in (-1, 1) for a risk difference')
    if (alternative == 'greater' && !(Delta > -margin)) {
      stop(name, " must exceed -margin when alternative is 'greater'")
    }
    if (alternative == 'less' && !(Delta < margin)) {
      stop(name, " must fall below margin when alternative is 'less'")
    }
    return(invisible(NULL))
  }
  if (effect == 'RD') {
    null <- 0
    if (abs(Delta) >= 1) stop(name, ' must lie in (-1, 1) for a risk difference')
  } else {
    null <- 1
    if (Delta <= 0) stop(name, ' must be positive for a risk ratio or an odds ratio')
  }
  if (Delta == null) stop(name, ' must differ from ', null)
  if (alternative == 'greater' && Delta < null) {
    stop(name, ' must exceed ', null, " when alternative is 'greater'")
  }
  if (alternative == 'less' && Delta > null) {
    stop(name, ' must fall below ', null, " when alternative is 'less'")
  }
  invisible(NULL)
}
