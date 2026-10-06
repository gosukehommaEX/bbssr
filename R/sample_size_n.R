#' Sample Sizes of the Two Groups for a Fixed Design
#'
#' Internal helper returning the sample sizes of a fixed design for given response
#' probabilities, either by a search over the exact power of the selected test or from
#' the normal approximation followed by the selected rounding rule. The arguments are
#' assumed to have been validated by the caller.
#'
#' Under \code{method = 'exact'} the search of \code{ss_exact_search} is used: it starts
#' from the normal approximation, lowers the size of group 2 one unit at a time as long
#' as the exact power attains the target power, or otherwise raises it one unit at a time
#' until the power attains it, and so returns a size that attains the target power while
#' the size one unit smaller does not (unless the returned size is 1). Group 1 receives
#' \code{ceiling(r N2)} patients. Under the two normal methods, the unrounded size
#' \code{n2} of group 2 is converted as follows, where \code{round} rounds halves up.
#' \describe{
#'   \item{group}{\code{N2 = ceiling(n2)} and \code{N1 = ceiling(r N2)}}
#'   \item{friede-kieser}{\code{N2 = ceiling(n2)} and \code{N1 = ceiling(r n2)}, as in
#'     Friede and Kieser (2004)}
#'   \item{total}{the total \code{N = ceiling((1 + r) n2)} is split into
#'     \code{N2 = floor(N / (1 + r))} and \code{N1 = N - N2}}
#'   \item{nearest}{\code{N2 = round(n2)} and \code{N1 = round(r n2)}, as in Farrington
#'     and Manning (1990)}
#' }
#' The ceilings of the normal methods ignore excesses below \code{1e-9}.
#'
#' @param p1 Response probability of group 1
#' @param p2 Response probability of group 2
#' @param r Allocation ratio to group 1
#' @param alpha Level of significance
#' @param tar.power Target power
#' @param Test Type of statistical test
#' @param alternative \code{'greater'}, \code{'less'} or \code{'two.sided'}
#' @param tsmethod \code{'minlike'}, \code{'central'} or \code{'blaker'}
#' @param n.grid Number of grid points over the nuisance parameter
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure
#' @param method \code{'exact'}, \code{'standard'} or \code{'null.variance'}
#' @param rounding \code{'group'}, \code{'friede-kieser'}, \code{'total'} or
#'   \code{'nearest'}
#' @param ref.pvalue Logical. Whether the maximum over the nuisance parameter of the
#'   unconditional tests is refined between the grid points
#' @param margin Non-inferiority margin on the scale of the risk difference, 0 for a test
#'   of superiority
#'
#' @return An integer vector with elements \code{N1} and \code{N2}
#'
#' @keywords internal
#' @noRd
sample_size_n <- function(p1, p2, r, alpha, tar.power, Test, alternative, tsmethod,
                          n.grid, bb.gamma, method, rounding, ref.pvalue, margin) {
  if (method != 'exact') {
    n2 <- ss_raw_n2(p1, p2, r, alpha, tar.power, alternative, method, margin)
    if (rounding == 'group') {
      N2 <- ceil_tol(n2)
      N1 <- ceiling(r * N2)
    } else if (rounding == 'friede-kieser') {
      N2 <- ceil_tol(n2)
      N1 <- ceil_tol(r * n2)
    } else if (rounding == 'nearest') {
      N2 <- floor(n2 + 0.5)
      N1 <- floor(r * n2 + 0.5)
    } else {
      N <- ceil_tol((1 + r) * n2)
      N2 <- floor(N / (1 + r))
      N1 <- N - N2
    }
    return(c(N1 = as.integer(max(1, N1)), N2 = as.integer(max(1, N2))))
  }
  N2 <- ss_exact_search(p1, p2, r, alpha, tar.power, Test, alternative, tsmethod, n.grid,
                        bb.gamma, ref.pvalue, margin)
  c(N1 = as.integer(ceiling(r * N2)), N2 = as.integer(N2))
}
