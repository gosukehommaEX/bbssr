#' Exact Power of a Stratified Test
#'
#' Calculates the exact power of a test of \code{\link{StratRR}} for a two-arm trial with a
#' binary endpoint and fixed group sizes in each stratum.
#'
#' @param p1 Response probabilities of group 1 in the strata: a vector with one element
#'   per stratum, or a matrix with one row per scenario and one column per stratum
#' @param p2 Response probabilities of group 2 in the strata, of the same form as
#'   \code{p1}
#' @inheritParams StratRR
#'
#' @return An object of class \code{bbssr_spower}, a data frame with one row per scenario
#'   containing the response probabilities \code{p1.S1, p1.S2, ...} and
#'   \code{p2.S1, p2.S2, ...} and the exact power \code{Power}.
#'
#' @details
#' The numbers of responders of the two groups in each stratum are independent binomial
#' variables. The power is the probability of the rejection region of
#' \code{\link{StratRR}}, computed exactly by enumerating the numbers of responders of the
#' strata and carrying the distribution of the partial sum of \code{S} over the strata.
#' The response probabilities may differ between the strata in any way, so the power is
#' available also when the odds ratio or the risk difference is not common to the strata.
#'
#' @examples
#' StratPower(p1 = c(0.5, 0.7), p2 = c(0.3, 0.5), N1 = c(40, 40), N2 = c(40, 40),
#'            alpha = 0.025, Test = 'CMH')
#'
#' # Several scenarios, and the stratified exact test
#' StratPower(p1 = rbind(c(0.5, 0.7), c(0.45, 0.65)), p2 = rbind(c(0.3, 0.5), c(0.3, 0.5)),
#'            N1 = c(40, 40), N2 = c(40, 40), alpha = 0.025, Test = 'Strat-Fisher')
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @seealso \code{\link{StratRR}}, \code{\link{StratSampleSize}}
#' @export
StratPower <- function(p1, p2, N1, N2, alpha, Test,
                       alternative = c('greater', 'less', 'two.sided'),
                       tsmethod = c('minlike', 'central', 'blaker')) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  a <- check_strat_args(N1, N2, alpha, Test, alternative, tsmethod)
  J <- length(a$N1)
  p1 <- check_strat_probs(p1, J, 'p1')
  p2 <- check_strat_probs(p2, J, 'p2')
  if (nrow(p1) != nrow(p2)) {
    if (nrow(p1) == 1) p1 <- p1[rep(1, nrow(p2)), , drop = FALSE]
    else if (nrow(p2) == 1) p2 <- p2[rep(1, nrow(p1)), , drop = FALSE]
    else stop('p1 and p2 must have the same number of scenarios')
  }
  Power <- strat_fixed_power(a$N1, a$N2, p1, p2, alpha, a$Test, a$alternative, a$tsmethod)
  out <- data.frame(p1, p2, Power = Power)
  names(out) <- c(paste0('p1.S', seq_len(J)), paste0('p2.S', seq_len(J)), 'Power')
  attr(out, 'N1') <- a$N1
  attr(out, 'N2') <- a$N2
  attr(out, 'alpha') <- alpha
  attr(out, 'Test') <- a$Test
  attr(out, 'alternative') <- a$alternative
  attr(out, 'tsmethod') <- a$tsmethod
  class(out) <- c('bbssr_spower', 'data.frame')
  out
}
