#' Exact Sample Size Search for Several Pairs of Response Probabilities
#'
#' Internal helper returning, for each pair of response probabilities, the size of
#' group 2 found by the search of \code{\link{BinarySampleSize}} with
#' \code{method = 'exact'}. Group 1 receives \code{ceiling(r N2)} patients.
#'
#' The search for each pair starts from the normal approximation and moves the size of
#' group 2 one unit at a time until the smallest size attaining the target power is
#' reached. The exact power at a candidate size is computed once for all pairs, by
#' \code{power_from_rr_multi}, and reused by every pair whose search visits that size. The
#' summation order is that of \code{power_from_rr}, so each pair receives the same sample
#' size as a separate search would give.
#'
#' @param p1 Response probabilities of group 1
#' @param p2 Response probabilities of group 2, of the same length as \code{p1}, with no
#'   element equal to the corresponding element of \code{p1}
#' @param r Allocation ratio to group 1
#' @param alpha Level of significance
#' @param tar.power Target power
#' @param Test Type of statistical test
#' @param alternative \code{'greater'}, \code{'less'} or \code{'two.sided'}
#' @param tsmethod \code{'minlike'} or \code{'central'}
#' @param n.grid Number of grid points over the nuisance parameter
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure
#'
#' @return An integer vector of sizes of group 2
#'
#' @keywords internal
#' @noRd
#' @import fpCompare
#' @importFrom stats qnorm dbinom
ss_exact_search <- function(p1, p2, r, alpha, tar.power, Test, alternative, tsmethod,
                            n.grid, bb.gamma) {
  store <- new.env(parent = emptyenv())
  # Exact power of every pair at a given size of group 2, computed once per size
  power_at <- function(N2) {
    key <- as.character(N2)
    hit <- store[[key]]
    if (!is.null(hit)) return(hit)
    N1 <- ceiling(r * N2)
    rr <- get_rr(N1, N2, alpha, Test, alternative, tsmethod, n.grid, bb.gamma)
    d1 <- outer(0:N1, p1, function(x, p) dbinom(x, N1, p))
    d2 <- outer(0:N2, p2, function(x, p) dbinom(x, N2, p))
    pw <- power_from_rr_multi(rr, d1, d2)
    assign(key, pw, envir = store)
    pw
  }
  # Step 0 (initial sample size from the normal approximation to the chi-squared test)
  alpha.eff <- if (alternative == 'two.sided') alpha / 2 else alpha
  p <- (r * p1 + p2) / (1 + r)
  init.N2 <- '*'(
    (1 + 1 / r) / ((p1 - p2) ^ 2),
    (qnorm(alpha.eff) * sqrt(p * (1 - p)) +
       qnorm(1 - tar.power) * sqrt((p1 * (1 - p1) / r + p2 * (1 - p2)) / (1 + 1 / r))) ^ 2
  )
  out <- integer(length(p1))
  for (k in seq_along(p1)) {
    # Step 1 (power calculation given the initial sample size)
    N2 <- max(1, ceiling(init.N2[k]))
    Power <- power_at(N2)[k]
    # Step 2 (sample size calculation via a grid search algorithm)
    if (Power %>=% tar.power) {
      while ((Power %>=% tar.power) && (N2 > 1)) {
        N2 <- N2 - 1
        Power <- power_at(N2)[k]
      }
      if (Power %<<% tar.power) N2 <- N2 + 1
    } else {
      while (Power %<<% tar.power) {
        N2 <- N2 + 1
        Power <- power_at(N2)[k]
      }
    }
    out[k] <- as.integer(N2)
  }
  out
}
