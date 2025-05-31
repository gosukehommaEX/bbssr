#' Clopper-Pearson Confidence Interval Bounds
#'
#' Computes the lower and upper bounds of Clopper-Pearson confidence intervals
#' for binomial proportions using the exact method based on beta distribution.
#'
#' @param x Number of successes (non-negative integer).
#' @param n Sample size (positive integer).
#' @param alpha Significance level for confidence interval (0 < alpha < 1).
#'
#' @return A numeric vector of length 2 containing the lower and upper bounds
#' of the confidence interval.
#'
#' @details
#' The Clopper-Pearson confidence interval is an exact confidence interval
#' for binomial proportions. It is based on the relationship between the
#' binomial distribution and the beta distribution. For x = 0 or x = n,
#' the bounds are set to the boundary values 0 or 1 respectively.
#'
#' The confidence interval has coverage probability of at least (1 - alpha),
#' making it a conservative but exact method.
#'
#' @examples
#' \dontrun{
#' # 95% confidence interval for 7 successes out of 20 trials
#' bounds <- ClopperPearsonCI(7, 20, 0.05)
#' print(bounds)
#'
#' # 99% confidence interval for edge cases
#' bounds_zero <- ClopperPearsonCI(0, 10, 0.01)
#' bounds_all <- ClopperPearsonCI(10, 10, 0.01)
#' }
#'
#' @references
#' Clopper, C. J., & Pearson, E. S. (1934). The use of confidence or fiducial
#' limits illustrated in the case of the binomial. Biometrika, 26(4), 404-413.
#'
#' @export
#' @import stats qbeta
ClopperPearsonCI <- function(x, n, alpha) {

  # Input validation
  if (n < 0 || x < 0 || x > n) {
    stop('Invalid input: n must be non-negative, x must be between 0 and n')
  }
  if (alpha <= 0 || alpha >= 1) {
    stop('alpha must be between 0 and 1')
  }

  # Handle edge case: n = 0
  if (n == 0) return(c(0, 1))

  # Lower bound
  if (x == 0) {
    lower <- 0
  } else {
    lower <- qbeta(alpha / 2, x, n - x + 1)
  }

  # Upper bound
  if (x == n) {
    upper <- 1
  } else {
    upper <- qbeta(1 - alpha / 2, x + 1, n - x)
  }

  return(c(lower, upper))
}
