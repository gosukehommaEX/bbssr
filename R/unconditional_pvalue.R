#' p-Values of an Unconditional Exact Test over the Whole Outcome Grid
#'
#' Internal helper shared by the Z-pooled and the Boschloo exact unconditional tests. Given
#' an ordering statistic evaluated at every cell of the outcome grid, the function returns
#' the null tail probability of that statistic maximized over the nuisance parameter, which
#' is the common response probability under the null hypothesis.
#'
#' Cells that share the same value of the ordering statistic form a tie group. The tail
#' event is defined as being at least as extreme as the observed cell, so every member of a
#' tie group receives the tail probability accumulated up to the last member of that group.
#' A plain cumulative sum over an arbitrary ordering within a tie group assigns different
#' p-values to cells that are statistically indistinguishable, which is why the tie
#' structure is computed explicitly.
#'
#' When \code{bb.gamma} is positive the Berger-Boos procedure is applied. The maximization
#' is then restricted to an exact \code{100 (1 - bb.gamma)} percent confidence interval for
#' the nuisance parameter, computed from the total number of responders of the observed
#' cell, and \code{bb.gamma} is added to the resulting tail probability. The confidence
#' bounds are appended to the grid so that the endpoints of the interval are always
#' evaluated.
#'
#' @param stat A numeric matrix of dimension \code{(N1 + 1)} by \code{(N2 + 1)} holding the
#'   ordering statistic
#' @param N1 Sample size for group 1
#' @param N2 Sample size for group 2
#' @param n.grid Number of grid points used to search over the nuisance parameter
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure. A value of 0
#'   disables the procedure
#' @param decreasing Logical. \code{TRUE} when a larger value of \code{stat} is more
#'   extreme, \code{FALSE} when a smaller value is more extreme
#'
#' @return A numeric matrix of p-values of the same dimension as \code{stat}
#'
#' @keywords internal
#' @noRd
#' @importFrom stats dbinom
unconditional_pvalue <- function(stat, N1, N2, n.grid, bb.gamma, decreasing) {
  N <- N1 + N2
  x1 <- c(row(stat)) - 1L
  x2 <- c(col(stat)) - 1L
  # Sort the cells from the most extreme to the least extreme
  ord <- order(c(stat), decreasing = decreasing)
  grp <- tie_groups(c(stat)[ord])
  # Grid of the nuisance parameter
  theta <- seq(0, 1, length.out = n.grid)
  if (bb.gamma > 0) {
    bnd <- cp_bounds(N, bb.gamma)
    theta <- sort(unique(c(theta, bnd[, 'lower'], bnd[, 'upper'])))
  }
  dbinom1 <- outer(0:N1, theta, function(x, t) dbinom(x, N1, t))
  dbinom2 <- outer(0:N2, theta, function(x, t) dbinom(x, N2, t))
  # Range of grid indices over which the tail probability is maximized, zero based
  if (bb.gamma > 0) {
    s.ord <- x1[ord] + x2[ord]
    g.lo <- findInterval(bnd[s.ord + 1L, 'lower'], theta) - 1L
    g.hi <- findInterval(bnd[s.ord + 1L, 'upper'], theta) - 1L
  } else {
    g.lo <- rep(0L, length(ord))
    g.hi <- rep(length(theta) - 1L, length(ord))
  }
  p.ord <- max_tail_prob(
    dbinom1, dbinom2,
    as.integer(x1[ord]), as.integer(x2[ord]),
    as.integer(grp[, 'last'] - 1L),
    as.integer(g.lo), as.integer(g.hi)
  )
  if (bb.gamma > 0) p.ord <- p.ord + bb.gamma
  p.val <- stat
  p.val[ord] <- pmin(1, p.ord)
  p.val
}
