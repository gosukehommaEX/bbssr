#' p-Values of the Exact Unconditional McNemar Test over the Whole Outcome Grid
#'
#' Internal helper returning the null tail probability of an ordering statistic,
#' maximized over the probability \code{psi = p10 + p01} of a discordant pair, at every
#' outcome of \code{N} pairs. Under the null hypothesis \code{p10 = p01 = psi / 2} the
#' number \code{D = a + b} of discordant pairs is binomial with \code{N} trials and
#' probability \code{psi}, and given \code{D} the number \code{a} is binomial with
#' \code{D} trials and probability 1/2, so that the tail probability of an outcome is a
#' polynomial in \code{psi} of degree \code{N}. This is the test of Suissa and Shuster
#' (1991) and the test Z of Berger and Sidik (2003), whose formula (8) takes the supremum
#' over the whole range of \code{psi}.
#'
#' Outcomes sharing the same value of the ordering statistic form a tie group and receive
#' the tail probability accumulated up to the last member of the group, as in
#' \code{unconditional_pvalue}. The maximum is taken over a grid of \code{n.grid} values of
#' \code{psi}. With \code{ref.pvalue = TRUE} the grid is extended by points equally spaced
#' on the arcsine square-root scale, as in \code{unconditional_pvalue}, and every local
#' maximum on the grid is refined by a golden-section search between its two neighbouring
#' grid points. With a positive \code{bb.gamma} the maximization is restricted to the
#' Clopper-Pearson interval of level \code{1 - bb.gamma} for \code{psi} computed from
#' \code{D}, whose end points are added to the grid, and \code{bb.gamma} is added to the
#' result, the confidence interval p-value of Berger and Sidik (2003, formula (13)).
#'
#' @param stat A numeric matrix of dimension \code{(N + 1)} by \code{(N + 1)} holding the
#'   ordering statistic, larger values being more extreme, with \code{NA} where
#'   \code{a + b > N}
#' @param N Number of pairs
#' @param n.grid Number of grid points over \code{psi}
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure. A value of 0
#'   disables the procedure
#' @param ref.pvalue Logical. Whether the maximum over \code{psi} is refined between the
#'   grid points
#'
#' @return A numeric matrix of p-values of the same dimension as \code{stat}, with
#'   \code{NA} where \code{a + b > N}
#'
#' @references
#' Berger RL, Sidik K (2003). Exact unconditional tests for a 2 x 2 matched-pairs design.
#' \emph{Statistical Methods in Medical Research}, 12(2), 91-108.
#'
#' Suissa S, Shuster JJ (1991). The 2 x 2 matched-pairs trial: exact unconditional design
#' and analysis. \emph{Biometrics}, 47(2), 361-372.
#'
#' @keywords internal
#' @noRd
#' @importFrom stats dbinom
paired_uncond_pvalue <- function(stat, N, n.grid, bb.gamma, ref.pvalue) {
  ok <- which(!is.na(stat))
  a <- row(stat)[ok] - 1L
  b <- col(stat)[ok] - 1L
  # Sort the outcomes from the most extreme to the least extreme
  ord <- order(stat[ok], decreasing = TRUE)
  grp <- tie_groups(stat[ok][ord])
  # Grid of psi
  theta <- seq(0, 1, length.out = n.grid)
  if (ref.pvalue) {
    phi <- seq(0, pi / 2, length.out = ceiling(4 * pi * sqrt(N)) + 1)
    theta <- sort(unique(c(theta, sin(phi) ^ 2)))
  }
  if (bb.gamma > 0) {
    bnd <- cp_bounds(N, bb.gamma)
    theta <- sort(unique(c(theta, bnd[, 'lower'], bnd[, 'upper'])))
  }
  dbinomD <- outer(0:N, theta, function(x, t) dbinom(x, N, t))
  D <- a[ord] + b[ord]
  w <- dbinom(a[ord], D, 0.5)
  # Range of grid indices over which the tail probability is maximized, zero based
  if (bb.gamma > 0) {
    g.lo <- findInterval(bnd[D + 1L, 'lower'], theta) - 1L
    g.hi <- findInterval(bnd[D + 1L, 'upper'], theta) - 1L
  } else {
    g.lo <- rep(0L, length(ord))
    g.hi <- rep(length(theta) - 1L, length(ord))
  }
  p.ord <- paired_tail_max(dbinomD, as.integer(D), w, as.integer(grp[, 'last'] - 1L),
                           as.integer(g.lo), as.integer(g.hi), theta, ref.pvalue)
  if (bb.gamma > 0) p.ord <- p.ord + bb.gamma
  p.val <- stat
  p.val[ok[ord]] <- pmin(1, p.ord)
  p.val
}
