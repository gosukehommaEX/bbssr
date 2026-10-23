#' Exact Power of a Fixed Stratified Design
#'
#' Internal helper returning the exact rejection probability of a stratified design with
#' the group sizes \code{N1} and \code{N2} for each row of \code{p1} and \code{p2}. With two
#' strata the stored critical values of \code{get_strat_crit} are used, through the
#' computation of \code{strat_bssr_power} for a design without an interim analysis; with
#' more strata \code{strat_power} enumerates the stratum totals and computes the critical
#' values as it goes.
#'
#' @param p1 Response probabilities of group 1, a matrix with one row per scenario and one
#'   column per stratum
#' @param p2 Response probabilities of group 2, of the same form
#' @inheritParams StratRR
#'
#' @return A numeric vector with one element per row of \code{p1}
#'
#' @keywords internal
#' @noRd
strat_fixed_power <- function(N1, N2, p1, p2, alpha, Test, alternative, tsmethod) {
  a <- check_strat_args(N1, N2, alpha, Test, alternative, tsmethod)
  J <- length(a$N1)
  p1 <- matrix(as.double(p1), ncol = J)
  p2 <- matrix(as.double(p2), ncol = J)
  if (J == 2) {
    crit <- get_strat_crit(a$N1, a$N2, alpha, a$Test, a$alternative, a$tsmethod)
    return(strat_bssr_power(c(0L, 0L), c(0L, 0L), list(crit), 0L, p1, p2))
  }
  strat_power(a$N1, a$N2, p1, p2, a$type, a$alt, a$ts, alpha)
}
