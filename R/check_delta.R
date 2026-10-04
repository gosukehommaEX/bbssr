#' Validate an Assumed Treatment Effect
#'
#' Internal helper checking that an assumed treatment effect is a single value on the
#' scale given by \code{effect}, that it differs from the value under the null hypothesis
#' and that its direction agrees with \code{alternative}.
#'
#' @param Delta Treatment effect
#' @param effect \code{'RD'}, \code{'RR'} or \code{'OR'}
#' @param alternative \code{'greater'}, \code{'less'} or \code{'two.sided'}
#' @param name Name of the argument, used in the error messages
#'
#' @return \code{NULL}, invisibly
#'
#' @keywords internal
#' @noRd
check_delta <- function(Delta, effect, alternative, name = 'Delta.A') {
  if (length(Delta) != 1 || is.na(Delta)) stop(name, ' must be a single value')
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
