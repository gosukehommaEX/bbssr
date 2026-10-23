#' Simulated Rejections of a Stratified Re-Estimation Design
#'
#' Internal helper simulating \code{nsim} trials of a stratified design with blinded sample
#' size re-estimation for one scenario. The interim numbers of responders of the two groups
#' of every stratum are drawn from binomial distributions, the final group sizes follow
#' from the pooled numbers of responders by \code{strat_reestimate}, the second-stage
#' numbers of responders are drawn given these sizes, and the final analysis is applied to
#' the combined data. The random numbers come from the generator of R, so the results
#' follow \code{set.seed()}.
#'
#' @param n1 Interim sizes of group 1 in the strata
#' @param n2 Interim sizes of group 2 in the strata
#' @param p1 Response probabilities of group 1 in the strata
#' @param p2 Response probabilities of group 2 in the strata
#' @param nsim Number of simulated trials
#' @param alpha Level of the final analysis
#' @param a List returned by \code{check_strat_args} for the final analysis
#' @param re.args List of the remaining arguments of \code{strat_reestimate}
#'
#' @return A list with the logical vector \code{reject} and the integer vector \code{N} of
#'   final total sample sizes of the simulated trials
#'
#' @keywords internal
#' @noRd
#' @importFrom stats rbinom
strat_simulate <- function(n1, n2, p1, p2, nsim, alpha, a, re.args) {
  J <- length(n1)
  draw <- function(size, prob) {
    x <- rbinom(length(size), size, prob)
    matrix(as.integer(x), ncol = J)
  }
  n1m <- matrix(n1, nsim, J, byrow = TRUE)
  n2m <- matrix(n2, nsim, J, byrow = TRUE)
  A <- draw(as.vector(n1m), rep(p1, each = nsim))
  B <- draw(as.vector(n2m), rep(p2, each = nsim))
  S <- A + B
  key <- apply(S, 1, paste, collapse = '_')
  first <- which(!duplicated(key))
  re <- do.call(strat_reestimate,
                c(list(s = S[first, , drop = FALSE], n1 = n1, n2 = n2,
                       alternative = a$alternative, tsmethod = a$tsmethod), re.args))
  idx <- match(key, key[first])
  N1 <- re$N1[idx, , drop = FALSE]
  N2 <- re$N2[idx, , drop = FALSE]
  X1 <- A + draw(as.vector(N1 - n1m), rep(p1, each = nsim))
  X2 <- B + draw(as.vector(N2 - n2m), rep(p2, each = nsim))
  pv <- strat_pvalue_cases(N1, N2, X1, X2, a$type, a$alt, a$ts)
  list(reject = pv %<<% alpha, N = as.integer(rowSums(N1) + rowSums(N2)))
}
