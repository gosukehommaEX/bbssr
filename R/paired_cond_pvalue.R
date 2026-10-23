#' p-Values of the Exact Conditional McNemar Test over the Whole Outcome Grid
#'
#' Internal helper returning the p-values of the exact conditional McNemar test, or its
#' mid-p version, at every outcome of \code{N} pairs. Given \code{D = a + b} discordant
#' pairs, the number \code{a} of pairs in which only treatment 1 responds is binomial with
#' \code{D} trials and probability 1/2 under the null hypothesis \code{p10 = p01}. For
#' \code{alternative = 'greater'} the p-value is \code{P(X >= a)} and its mid-p version
#' \code{P(X > a) + f(a) / 2}, with \code{X} binomial with \code{D} trials and probability
#' 1/2 and \code{f} its probability function. The two-sided p-value is twice the smaller
#' tail probability, at most 1, formula (4) of Fagerland et al. (2013). The two-sided mid-p
#' value is the two-sided p-value minus \code{f(a)}, or \code{1 - f(a) / 2} when
#' \code{a = b}, formulae (5) and (6) of that article. The conditional distribution is
#' symmetric, so the conventions \code{'minlike'}, \code{'central'} and \code{'blaker'}
#' give the same two-sided p-value, and the mid-p value of Fagerland et al. (2013) is
#' that of \code{'minlike'} and \code{'blaker'} with the tied outcomes contributing half
#' of their probability. The alternative \code{'less'} is obtained by exchanging \code{a}
#' and \code{b}.
#'
#' @param N Number of pairs
#' @param alternative \code{'greater'} or \code{'two.sided'}
#' @param midp Logical. Whether the mid-p value is returned
#'
#' @return A numeric matrix of dimension \code{(N + 1)} by \code{(N + 1)} whose element in
#'   row \code{a + 1} and column \code{b + 1} is the p-value of the outcome
#'   \code{(a, b)}, and \code{NA} where \code{a + b > N}
#'
#' @references
#' Fagerland MW, Lydersen S, Laake P (2013). The McNemar test for binary matched-pairs
#' data: mid-p and asymptotic are better than exact conditional. \emph{BMC Medical
#' Research Methodology}, 13, 91.
#'
#' @keywords internal
#' @noRd
#' @importFrom stats dbinom pbinom
paired_cond_pvalue <- function(N, alternative, midp) {
  P <- matrix(NA_real_, N + 1, N + 1)
  ok <- outer(0:N, 0:N, '+') <= N
  a <- row(P)[ok] - 1
  b <- col(P)[ok] - 1
  D <- a + b
  f <- dbinom(a, D, 0.5)
  # Upper and lower tail probabilities P(X >= a) and P(X <= a)
  up <- pbinom(a - 1, D, 0.5, lower.tail = FALSE)
  lo <- pbinom(a, D, 0.5)
  if (alternative == 'greater') {
    p <- if (midp) up - f / 2 else up
  } else {
    p <- pmin(1, 2 * pmin(lo, up))
    if (midp) p <- ifelse(a == b, 1 - f / 2, 2 * pmin(lo, up) - f)
  }
  P[ok] <- p
  P
}
