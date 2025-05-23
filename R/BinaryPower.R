#' Power Calculation for Two-Arm Trials with Binary Endpoints
#'
#' Calculates power for two-arm trials with binary endpoints using exact statistical tests.
#' The function supports five different one-sided tests and can handle vectors of probabilities.
#'
#' @param p1 True probability of responders for group 1 (can be a vector with different values)
#' @param p2 True probability of responders for group 2 (can be a vector with different values)
#' @param N1 Sample size for group 1
#' @param N2 Sample size for group 2
#' @param alpha One-sided level of significance
#' @param Test Type of statistical test. Options: 'Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', or 'Boschloo'
#'
#' @return A numeric value or vector of power values. If vectors are provided for p1 and p2,
#' a vector of powers corresponding to each combination will be returned.
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
#' The power calculation is based on the exact distribution of the test statistic
#' under the specified alternative hypothesis.
#'
#' @examples
#' \dontrun{
#' # Single power calculation
#' power1 <- BinaryPower(p1 = 0.5, p2 = 0.2, N1 = 10, N2 = 40,
#'                      alpha = 0.025, Test = 'Chisq')
#' print(power1)
#'
#' # Multiple power calculations
#' p1_vec <- c(0.5, 0.6, 0.7, 0.8)
#' p2_vec <- c(0.2, 0.2, 0.2, 0.2)
#' powers <- BinaryPower(p1 = p1_vec, p2 = p2_vec, N1 = 10, N2 = 40,
#'                      alpha = 0.025, Test = 'Fisher')
#' print(powers)
#' }
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @export
#' @import fpCompare
#' @import stats
BinaryPower <- function(p1, p2, N1, N2, alpha, Test) {
  # Check that p1 and p2 are the same length
  if(length(p1) %!=% length(p2)) stop('p1 and p2 should be the same length')
  # Set rejection region
  RR <- BinaryRR(N1, N2, alpha, Test)
  # Return power
  sapply(seq(length(p1)), function(i) sum(dbinom(0:N1, N1, p1[i]) * pbinom(rowSums(RR) - 1, N2, p2[i])))
}
