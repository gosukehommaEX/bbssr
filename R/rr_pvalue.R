#' p-Values of a Test over the Whole Outcome Grid
#'
#' Internal helper returning the p-value of the selected test at every cell of the
#' \code{(N1 + 1)} by \code{(N2 + 1)} outcome grid. The arguments are assumed to have been
#' validated by \code{check_rr_args}.
#'
#' A lower-tail alternative is obtained from the upper-tail version with the two groups
#' exchanged. Every test considered treats the groups symmetrically apart from the
#' direction of the alternative, so the p-value of the outcome \code{(x1, x2)} against
#' \code{p1 < p2} equals the p-value of \code{(x2, x1)} against \code{p1 > p2} in a trial
#' with the sample sizes exchanged.
#'
#' @param N1 Sample size for group 1
#' @param N2 Sample size for group 2
#' @param Test Full name of the test
#' @param alternative \code{'greater'}, \code{'less'} or \code{'two.sided'}
#' @param tsmethod \code{'minlike'} or \code{'central'}
#' @param n.grid Number of grid points over the nuisance parameter
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure
#' @param ref.pvalue Logical. Whether the maximum over the nuisance parameter of the
#'   unconditional tests is refined between the grid points
#' @param margin Non-inferiority margin on the scale of the risk difference, 0 for a test
#'   of superiority
#'
#' @return A numeric matrix of dimension \code{(N1 + 1)} by \code{(N2 + 1)}
#'
#' @keywords internal
#' @noRd
#' @importFrom stats pnorm
rr_pvalue <- function(N1, N2, Test, alternative, tsmethod, n.grid, bb.gamma, ref.pvalue,
                      margin) {
  if (alternative == 'less') {
    return(t(rr_pvalue(N2, N1, Test, 'greater', tsmethod, n.grid, bb.gamma, ref.pvalue,
                       margin)))
  }
  if (Test == 'Chisq') {
    Z <- zstat(N1, N2)
    if (alternative == 'greater') {
      pnorm(Z, lower.tail = FALSE)
    } else {
      pmin(2 * pnorm(abs(Z), lower.tail = FALSE), 1)
    }
  } else if (Test == 'Fisher') {
    fisher_pvalue(N1, N2, alternative, tsmethod, midp = FALSE)
  } else if (Test == 'Fisher-midP') {
    fisher_pvalue(N1, N2, alternative, tsmethod, midp = TRUE)
  } else if (Test %in% c('Blackwelder', 'Farrington-Manning')) {
    se <- if (Test == 'Blackwelder') 'unpooled' else 'restricted'
    Z <- zstat_margin(N1, N2, margin, se)
    if (alternative == 'greater') {
      pnorm(Z, lower.tail = FALSE)
    } else {
      pmin(2 * pnorm(abs(Z), lower.tail = FALSE), 1)
    }
  } else if (Test == 'Z-pool') {
    Z <- zstat(N1, N2)
    stat <- if (alternative == 'greater') Z else abs(Z)
    unconditional_pvalue(stat, N1, N2, n.grid, bb.gamma, decreasing = TRUE, ref.pvalue)
  } else {
    stat <- fisher_pvalue(N1, N2, alternative, tsmethod, midp = FALSE)
    unconditional_pvalue(stat, N1, N2, n.grid, bb.gamma, decreasing = FALSE, ref.pvalue)
  }
}
