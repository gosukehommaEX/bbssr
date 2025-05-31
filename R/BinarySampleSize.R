#' Sample Size Calculation for Two-Arm Trials with Binary Endpoints
#'
#' Calculates the required sample size for two-arm trials with binary endpoints using
#' various exact statistical tests. The function supports five different one-sided tests.
#' This optimized version provides consistent results and improved computational efficiency.
#'
#' @param p1 True probability of responders for group 1.
#' @param p2 True probability of responders for group 2.
#' @param r Allocation ratio to group 1 (i.e., allocation ratio of group 1:group 2 = r:1, r > 0).
#' @param alpha One-sided level of significance.
#' @param tar.power Target power.
#' @param Test Type of statistical test. Options: 'Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', or 'Boschloo'.
#' @param BB Logical. If TRUE, apply Berger-Boos approach for 'Z-pool' and 'Boschloo' tests (default: FALSE).
#' @param gamma Berger-Boos parameter for confidence interval adjustment (default: 0.0001).
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
#' The optimized calculation uses a three-step approach:
#' \enumerate{
#'   \item Calculate conservative initial sample size (intentionally smaller than required)
#'   \item Use binary search for efficient range exploration when power gap is large
#'   \item Fine-tune with sequential search for precise final result
#' }
#'
#' This approach ensures consistent results by always searching in the increasing direction
#' and provides significant speed improvements through intelligent search algorithms.
#'
#' When BB = TRUE, the Berger-Boos approach based on Clopper-Pearson confidence intervals
#' is applied to 'Z-pool' and 'Boschloo' tests for improved computational efficiency while
#' maintaining statistical validity.
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
#'
#' # Sample size for Boschloo test with Berger-Boos approach
#' result3 <- BinarySampleSize(p1 = 0.6, p2 = 0.3, r = 1, alpha = 0.025,
#'                            tar.power = 0.8, Test = 'Boschloo',
#'                            BB = TRUE, gamma = 0.0001)
#' print(result3)
#' }
#'
#' @export
#' @import fpCompare
#' @import stats qnorm
BinarySampleSize <- function(p1, p2, r, alpha, tar.power, Test, BB = FALSE, gamma = 0.0001) {

  # Input validation
  if (p1 < 0 || p1 > 1) stop('p1 must be between 0 and 1')
  if (p2 < 0 || p2 > 1) stop('p2 must be between 0 and 1')
  if (r <= 0) stop('r must be positive')
  if (alpha <= 0 || alpha >= 1) stop('alpha must be between 0 and 1')
  if (tar.power <= 0 || tar.power >= 1) stop('tar.power must be between 0 and 1')
  if (p1 <= p2) stop('p1 must be greater than p2 for one-sided test')

  # Step 0: Calculate conservative initial sample size
  # Use a multiplier < 1 to ensure we start below the required sample size
  p <- (r * p1 + p2) / (1 + r)
  conservative_multiplier <- 0.85  # Start with 85% of normal approximation

  init_N2 <- '*'(
    conservative_multiplier * (1 + 1 / r) / ((p1 - p2) ^ 2),
    (qnorm(alpha) * sqrt(p * (1 - p)) + qnorm(1 - tar.power) * sqrt((p1 * (1 - p1) / r + p2 * (1 - p2)) / (1 + 1 / r))) ^ 2
  )

  # Ensure minimum sample size of 2 for each group
  N2 <- max(2, floor(init_N2))
  N1 <- max(2, ceiling(r * N2))

  # Step 1: Initial power calculation
  Power <- BinaryPower(p1, p2, N1, N2, alpha, Test, BB = BB, gamma = gamma)

  # Step 2: Efficient search algorithm
  if (Power %>=% tar.power) {
    # Rare case: initial estimate was too high, decrease until below target
    while (Power %>=% tar.power && N2 > 2) {
      N2 <- N2 - 1
      N1 <- max(2, ceiling(r * N2))
      Power <- BinaryPower(p1, p2, N1, N2, alpha, Test, BB = BB, gamma = gamma)
    }
    # Now increase to meet target (ensures consistency)
    while (Power %<<% tar.power) {
      N2 <- N2 + 1
      N1 <- ceiling(r * N2)
      Power <- BinaryPower(p1, p2, N1, N2, alpha, Test, BB = BB, gamma = gamma)
    }
  } else {
    # Expected case: power is below target, need to increase

    # Phase 2a: Binary search for large gaps (speed optimization)
    if (Power < tar.power - 0.1) {  # Only use binary search if gap is substantial

      lower_N2 <- N2
      upper_N2 <- N2

      # Find upper bound using exponential search
      while (Power %<<% tar.power) {
        lower_N2 <- upper_N2
        upper_N2 <- upper_N2 * 2  # Exponential increase
        N1 <- ceiling(r * upper_N2)
        Power <- BinaryPower(p1, p2, N1, upper_N2, alpha, Test, BB = BB, gamma = gamma)
      }

      # Binary search between lower_N2 and upper_N2
      while (upper_N2 - lower_N2 > 5) {  # Switch to linear search when close
        mid_N2 <- floor((lower_N2 + upper_N2) / 2)
        N1 <- ceiling(r * mid_N2)
        mid_Power <- BinaryPower(p1, p2, N1, mid_N2, alpha, Test, BB = BB, gamma = gamma)

        if (mid_Power %<<% tar.power) {
          lower_N2 <- mid_N2
        } else {
          upper_N2 <- mid_N2
        }
      }

      # Set N2 to lower bound for final linear search
      N2 <- lower_N2
      N1 <- ceiling(r * N2)
      Power <- BinaryPower(p1, p2, N1, N2, alpha, Test, BB = BB, gamma = gamma)
    }

    # Phase 2b: Linear search for final precision
    while (Power %<<% tar.power) {
      N2 <- N2 + 1
      N1 <- ceiling(r * N2)
      Power <- BinaryPower(p1, p2, N1, N2, alpha, Test, BB = BB, gamma = gamma)
    }
  }

  # Step 3: Finalize results
  N <- N1 + N2

  # Final power calculation to ensure accuracy
  Power <- BinaryPower(p1, p2, N1, N2, alpha, Test, BB = BB, gamma = gamma)

  # Return result
  result <- data.frame(p1, p2, r, alpha, tar.power, Test, BB, gamma, Power, N1, N2, N)
  return(result)
}
