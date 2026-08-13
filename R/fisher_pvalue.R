#' Conditional Exact p-Values over the Whole Outcome Grid
#'
#' Internal helper returning the Fisher exact or Fisher mid-p p-value at every cell of the
#' \code{(N1 + 1)} by \code{(N2 + 1)} outcome grid. The p-values are obtained from the
#' hypergeometric distribution of the group 1 responder count conditional on the total
#' number of responders.
#'
#' For a two-sided alternative two conventions are available. The \code{minlike}
#' convention sums the null probabilities of all tables whose probability does not exceed
#' the probability of the observed table, which is the convention of
#' \code{stats::fisher.test}. The \code{central} convention doubles the smaller of the two
#' one-sided tail probabilities and truncates the result at 1. Under the mid-p correction
#' the observed table contributes half of its probability instead of its full probability.
#'
#' @param N1 Sample size for group 1
#' @param N2 Sample size for group 2
#' @param alternative Either \code{'greater'} or \code{'two.sided'}
#' @param tsmethod Either \code{'minlike'} or \code{'central'}, used only when
#'   \code{alternative} is \code{'two.sided'}
#' @param midp Logical. If \code{TRUE}, the mid-p correction is applied
#'
#' @return A numeric matrix of dimension \code{(N1 + 1)} by \code{(N2 + 1)}
#'
#' @keywords internal
#' @noRd
#' @importFrom stats dhyper
fisher_pvalue <- function(N1, N2, alternative, tsmethod, midp) {
  p.val <- matrix(NA_real_, nrow = N1 + 1, ncol = N2 + 1)
  for (s in 0:(N1 + N2)) {
    k <- max(0L, s - N2):min(N1, s)
    d <- dhyper(k, N1, N2, s)
    if (alternative == 'greater') {
      # One-sided upper tail probability P(X1 >= k)
      upper <- rev(cumsum(rev(d)))
      p.k <- if (midp) upper - 0.5 * d else upper
    } else if (tsmethod == 'central') {
      upper <- rev(cumsum(rev(d)))
      lower <- cumsum(d)
      if (midp) {
        upper <- upper - 0.5 * d
        lower <- lower - 0.5 * d
      }
      p.k <- 2 * pmin(lower, upper)
    } else {
      # Sum the probabilities of all tables that are no more likely than the observed one
      ord <- order(d)
      d.ord <- d[ord]
      cum.d <- cumsum(d.ord)
      grp <- tie_groups(d.ord)
      cum.d0 <- c(0, cum.d)
      p.ord <- if (midp) {
        0.5 * (cum.d[grp[, 'last']] + cum.d0[grp[, 'first']])
      } else {
        cum.d[grp[, 'last']]
      }
      p.k <- numeric(length(d))
      p.k[ord] <- p.ord
    }
    p.val[cbind(k + 1L, s - k + 1L)] <- pmin(1, p.k)
  }
  p.val
}
