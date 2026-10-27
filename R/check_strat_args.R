#' Validate the Arguments of a Stratified Rejection Region
#'
#' Internal helper checking the group sizes of the strata, the level, the test and the
#' alternative of a stratified design, and returning them together with the integer codes
#' of the test, the alternative and the two-sided method used by the compiled code.
#'
#' @param N1 Sizes of group 1 in the strata, one element per stratum
#' @param N2 Sizes of group 2 in the strata, of the same length
#' @param alpha Level of significance
#' @param Test \code{'CMH'}, \code{'CMH-cc'}, \code{'Cochran'}, \code{'Cochran-cc'},
#'   \code{'Strat-Fisher'}, \code{'Unstrat-Chisq'}, \code{'Unstrat-Chisq-cc'} or
#'   \code{'Unstrat-Fisher'}, partially matched
#' @param alternative \code{'greater'}, \code{'less'} or \code{'two.sided'}
#' @param tsmethod \code{'minlike'}, \code{'central'} or \code{'blaker'}
#'
#' @return A list with the integer vectors \code{N1} and \code{N2}, the full name
#'   \code{Test}, the codes \code{type}, \code{alt} and \code{ts}, \code{cc}, whether
#'   the test has a continuity correction, and \code{unstrat}, whether it ignores the
#'   strata. The two-sided method is set to \code{'minlike'} unless the test is
#'   \code{'Strat-Fisher'} or \code{'Unstrat-Fisher'} with a two-sided alternative, where
#'   it matters, so that equal regions share one code
#'
#' @keywords internal
#' @noRd
check_strat_args <- function(N1, N2, alpha, Test, alternative, tsmethod) {
  tests <- c('CMH', 'CMH-cc', 'Cochran', 'Cochran-cc', 'Strat-Fisher', 'Unstrat-Chisq',
             'Unstrat-Chisq-cc', 'Unstrat-Fisher')
  Test <- match.arg(Test, tests)
  if (!is.numeric(N1) || !is.numeric(N2) || length(N1) != length(N2) || length(N1) < 2 ||
      any(!is.finite(N1)) || any(!is.finite(N2)) || any(N1 < 1) || any(N2 < 1) ||
      any(N1 != round(N1)) || any(N2 != round(N2)) ||
      max(N1, N2) > .Machine$integer.max) {
    stop('N1 and N2 must be vectors of positive integers of the same length, at least ',
         '2, with one element per stratum')
  }
  if (length(alpha) != 1 || !is.numeric(alpha) || is.na(alpha) || alpha <= 0 ||
      alpha >= 1) {
    stop('alpha must be a single value in (0, 1)')
  }
  alternative <- match.arg(alternative, c('greater', 'less', 'two.sided'))
  tsmethod <- match.arg(tsmethod, c('minlike', 'central', 'blaker'))
  if (!(Test %in% c('Strat-Fisher', 'Unstrat-Fisher')) || alternative != 'two.sided') {
    tsmethod <- 'minlike'
  }
  list(N1 = as.integer(N1), N2 = as.integer(N2), Test = Test, alternative = alternative,
       tsmethod = tsmethod,
       type = match(Test, tests) - 1L,
       alt = match(alternative, c('greater', 'less', 'two.sided')) - 1L,
       ts = match(tsmethod, c('minlike', 'central', 'blaker')) - 1L,
       cc = Test %in% c('CMH-cc', 'Cochran-cc', 'Unstrat-Chisq-cc'),
       unstrat = Test %in% c('Unstrat-Chisq', 'Unstrat-Chisq-cc', 'Unstrat-Fisher'))
}
