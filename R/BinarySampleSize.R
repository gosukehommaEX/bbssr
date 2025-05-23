#' Sample Size Calculation for Two-Arm Trials with Binary Endpoints
#'
#' Calculates the required sample size for two-arm trials with binary endpoints using
#' various exact statistical tests. The function supports five different one-sided tests.
#'
#' @param p1 True probability of responders for group 1
#' @param p2 True probability of responders for group 2
#' @param r Allocation ratio to group 1 (i.e., allocation ratio of group 1:group 2 = r:1, r > 0)
#' @param alpha One-sided level of significance
#' @param tar.power Target power
#' @param Test Type of statistical test. Options: 'Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', or 'Boschloo'
#'
#' @return A data frame containing:
#' \describe{
#'   \item{p1}{True probability of responders for group 1}
#'   \item{p2}{True probability of responders for group 2}
#'   \item{r}{Allocation ratio to group 1}
#'   \item{alpha}{One-sided level of significance}
#'   \item{tar.power}{Target power}
#'   \item{Test}{Name of the statistical test}
#'   \item{Power}{Calculated power}
#'   \item{N1}{Required sample size of group 1}
#'   \item{N2}{Required sample size of group 2}
#'   \item{N}{Total required sample size}
#' }
#'
#' @details
#' The function supports the following five one-sided tests:
#' \itemize{
#'   \item The one-sided Pearson chi-squared test (Chisq)
#'   \item The Fisher exact test (Fisher)
#'   \item The Fisher mid-p test (Fisher-midP)
#'   \item The Z-pooled exact unconditional test (Z-pool)
#'   \item The Boschloo exact unconditional test (Boschloo)
#' }
#'
#' The calculation uses a three-step approach:
#' \enumerate{
#'   \item Calculate initial sample size using normal approximation for chi-squared test
#'   \item Perform power calculation with the initial sample size
#'   \item Use grid search algorithm to find the optimal sample size
#' }
#'
#' @examples
#' \dontrun{
#' # Sample size for Pearson chi-squared test
#' result1 <- BinarySampleSize(p1 = 0.4, p2 = 0.2, r = 2, alpha = 0.025,
#'                            tar.power = 0.8, Test = 'Chisq')
#' print(result1)
#'
#' # Sample size for Fisher exact test
#' result2 <- BinarySampleSize(p1 = 0.5, p2 = 0.2, r = 3, alpha = 0.025,
#'                            tar.power = 0.9, Test = 'Fisher')
#' print(result2)
#' }
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @export
#' @import fpCompare
#' @import stats
BinarySampleSize <- function(p1, p2, r, alpha, tar.power, Test) {
  # Step 0 (calculate the required sample size for the one-sided Pearson chi-squared test based on the normal approximation)
  p <- (r * p1 + p2) / (1 + r)
  init_N2 <- '*'(
    (1 + 1 / r) / ((p1 - p2) ^ 2),
    (qnorm(alpha) * sqrt(p * (1 - p)) + qnorm(1 - tar.power) * sqrt((p1 * (1 - p1) / r + p2 * (1 - p2)) / (1 + 1 / r))) ^ 2
  )
  # Step 1 (power calculation given initial sample size)
  N2 <- ceiling(init_N2)
  N1 <- ceiling(r * N2)
  Power <- BinaryPower(p1, p2, N1, N2, alpha, Test)
  # Step 2 (sample size calculation via a grid search algorithm)
  if(Power %>=% tar.power) {
    while(Power %>=% tar.power) {
      N2 <- N2 - 1
      N1 <- ceiling(r * N2)
      Power <- BinaryPower(p1, p2, N1, N2, alpha, Test)
    }
    N2 <- N2 + 1
  } else {
    while(Power %<<% tar.power) {
      N2 <- N2 + 1
      N1 <- ceiling(r * N2)
      Power <- BinaryPower(p1, p2, N1, N2, alpha, Test)
    }
  }
  # Step 3 (determine the final sample size)
  N1 <- ceiling(r * N2)
  N <- N1 + N2
  Power <- BinaryPower(p1, p2, N1, N2, alpha, Test)
  # Return result
  result <- data.frame(p1, p2, r, alpha, tar.power, Test, Power, N1, N2, N)
  return(result)
}
