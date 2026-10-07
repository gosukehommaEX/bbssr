#' Re-Estimated Sample Sizes for Recovered Response Probabilities
#'
#' Internal helper returning, for each pair of recovered response probabilities, the
#' sample sizes that a fixed design would require, together with the unrounded total
#' when the normal approximation is used. Each distinct pair is evaluated once.
#'
#' A pair of equal probabilities admits no sample size. The risk ratio produces one when
#' no interim patient responds, and the odds ratio when no interim patient or every
#' interim patient responds. With a non-inferiority margin, a pair on the boundary of the
#' null hypothesis admits no sample size instead. The planned sizes
#' \code{N1.plan} and \code{N2.plan} are returned for such a pair.
#'
#' @param hat.p1 Recovered response probabilities of group 1. They must lie in the unit
#'   interval when \code{method = 'exact'}, and are used before truncation otherwise
#' @param hat.p2 Recovered response probabilities of group 2, as \code{hat.p1}
#' @param r Allocation ratio to group 1
#' @param alpha Level of significance used for the re-estimation
#' @param tar.power Target power
#' @param Test Test whose exact power is used when \code{method = 'exact'}
#' @param alternative \code{'greater'}, \code{'less'} or \code{'two.sided'}
#' @param tsmethod \code{'minlike'}, \code{'central'} or \code{'blaker'}
#' @param n.grid Number of grid points over the nuisance parameter
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure
#' @param method \code{'exact'}, \code{'standard'} or \code{'null.variance'}
#' @param rounding \code{'group'}, \code{'friede-kieser'} or \code{'total'}
#' @param N1.plan Planned sample size of group 1, or \code{NULL}
#' @param N2.plan Planned sample size of group 2, or \code{NULL}
#' @param ref.pvalue Logical. Whether the maximum over the nuisance parameter of the
#'   unconditional tests is refined between the grid points
#' @param margin Non-inferiority margin on the scale of the risk difference, 0 for a test
#'   of superiority
#' @param search Search of \code{ss_exact_search} used under \code{method = 'exact'}
#' @param search.limit Limit of \code{search = 'stable'}, see \code{ss_exact_search}
#'
#' @return A data frame with columns \code{n.raw} (unrounded total, \code{NA} under
#'   \code{method = 'exact'}), \code{N1.re}, \code{N2.re} and \code{N2.limit}, the
#'   limit of the size of group 2 examined by \code{search = 'stable'} (otherwise
#'   \code{NA})
#'
#' @keywords internal
#' @noRd
#' @import fpCompare
reestimate <- function(hat.p1, hat.p2, r, alpha, tar.power, Test, alternative, tsmethod,
                       n.grid, bb.gamma, method, rounding, N1.plan, N2.plan, ref.pvalue,
                       margin, search = 'crossing', search.limit = c(2, 50)) {
  degenerate <- if (margin == 0) {
    hat.p1 %==% hat.p2
  } else if (alternative == 'less') {
    (hat.p2 - hat.p1 + margin) %==% 0
  } else {
    (hat.p1 - hat.p2 + margin) %==% 0
  }
  if (any(degenerate) && (is.null(N1.plan) || is.null(N2.plan))) {
    stop(if (margin == 0) 'the recovered response probabilities coincide' else
           'the recovered response probabilities lie on the null boundary',
         ', so no sample size can be re-estimated; supply the planned sample sizes N1 ',
         'and N2')
  }
  n <- length(hat.p1)
  n.raw <- rep(NA_real_, n)
  N1.re <- integer(n)
  N2.re <- integer(n)
  N2.limit <- rep(NA_integer_, n)
  L.u <- integer(0)
  if (method != 'exact') {
    n.raw <- '*'(1 + r, ss_raw_n2(hat.p1, hat.p2, r, alpha, tar.power, alternative, method,
                                  margin))
  }
  key <- paste(hat.p1, hat.p2, sep = '_')
  first <- which(!duplicated(key) & !degenerate)
  if (method == 'exact') {
    # One search over all distinct pairs, sharing the power at every candidate size
    s <- ss_exact_search(hat.p1[first], hat.p2[first], r, alpha, tar.power, Test,
                         alternative, tsmethod, n.grid, bb.gamma, ref.pvalue, margin,
                         search, search.limit)
    N2.u <- s$N2
    L.u <- s$limit
    N1.u <- as.integer(ceiling(r * N2.u))
  } else if (length(first) > 0) {
    n.u <- vapply(first, function(u) {
      sample_size_n(hat.p1[u], hat.p2[u], r, alpha, tar.power, Test, alternative,
                    tsmethod, n.grid, bb.gamma, method, rounding, ref.pvalue, margin)
    }, integer(2))
    N1.u <- n.u['N1', ]
    N2.u <- n.u['N2', ]
  } else {
    N1.u <- integer(0)
    N2.u <- integer(0)
  }
  idx <- match(key, key[first])
  ok <- !is.na(idx)
  N1.re[ok] <- N1.u[idx[ok]]
  N2.re[ok] <- N2.u[idx[ok]]
  if (length(L.u) > 0) N2.limit[ok] <- L.u[idx[ok]]
  if (any(degenerate)) {
    N1.re[degenerate] <- as.integer(N1.plan)
    N2.re[degenerate] <- as.integer(N2.plan)
    if (method != 'exact') n.raw[degenerate] <- N1.plan + N2.plan
  }
  data.frame(n.raw = n.raw, N1.re = N1.re, N2.re = N2.re, N2.limit = N2.limit)
}
