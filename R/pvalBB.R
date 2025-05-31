#' Optimized P-value Computation with Berger-Boos Approach
#'
#' Computes p-values using the Berger-Boos approach for exact unconditional tests.
#' This function provides computational efficiency while maintaining statistical
#' accuracy by restricting the nuisance parameter space. It works for both Z-pooled
#' and Boschloo tests based on the pre-sorted input order.
#'
#' @param x1 Vector of success counts for group 1 (pre-sorted by test statistic).
#' @param x2 Vector of success counts for group 2 (pre-sorted by test statistic).
#' @param N1 Sample size for group 1.
#' @param N2 Sample size for group 2.
#' @param gamma Berger-Boos parameter for confidence interval adjustment.
#'
#' @return A numeric vector of p-values corresponding to each (x1,x2) pair
#'
#' @details
#' This function implements the Berger-Boos approach for computing exact
#' unconditional p-values efficiently. It uses Clopper-Pearson confidence
#' intervals to restrict the theta parameter space, then computes the
#' maximum p-value over this restricted range. The final p-value includes
#' the gamma adjustment as per the Berger-Boos methodology.
#'
#' The function assumes that the input vectors x1 and x2 are pre-sorted
#' according to the test statistic (Z-statistic for Z-pool test in descending
#' order, or Fisher p-values for Boschloo test in ascending order). This
#' allows the function to work correctly for both test types without needing
#' to know which specific test is being used.
#'
#' The function uses vectorized operations and optimized matrix computations
#' to achieve high performance while maintaining consistency with non-BB
#' approaches by using the same theta evaluation grid.
#'
#' @examples
#' \dontrun{
#' # For Z-pool test: x1, x2 should be sorted by Z-statistic (descending)
#' x1 <- c(10, 9, 8)  # Pre-sorted by Z-statistic
#' x2 <- c(4, 3, 2)   # Pre-sorted by Z-statistic
#' p.vals.BB <- pvalBB(x1, x2, 20, 15, 0.001)
#' print(p.vals.BB)
#'
#' # For Boschloo test: x1, x2 should be sorted by Fisher p-value (ascending)
#' x1 <- c(8, 9, 10)  # Pre-sorted by Fisher p-value
#' x2 <- c(2, 3, 4)   # Pre-sorted by Fisher p-value
#' p.vals.BB <- pvalBB(x1, x2, 20, 15, 0.001)
#' print(p.vals.BB)
#' }
#'
#' @references
#' Berger, R. L., & Boos, D. D. (1994). P values maximized over a confidence
#' set for the nuisance parameter. Journal of the American Statistical
#' Association, 89(427), 1012-1016.
#'
#' @export
#' @import stats dbinom
pvalBB <- function(x1, x2, N1, N2, gamma) {
  # Input validation
  if (length(x1) != length(x2)) {
    stop('x1 and x2 must have the same length')
  }
  if (any(x1 < 0) || any(x2 < 0)) {
    stop('All values in x1 and x2 must be non-negative')
  }
  if (any(x1 > N1) || any(x2 > N2)) {
    stop('Values in x1 and x2 cannot exceed their respective sample sizes')
  }
  if (gamma <= 0 || gamma >= 1) {
    stop('gamma must be between 0 and 1')
  }

  n_points <- length(x1)

  # Use consistent theta grid with non-BB approach
  theta_seq <- seq(0, 1, length.out = 100)

  # Precompute unique binomial probabilities for efficiency
  uniq_i <- sort(unique(x1))
  uniq_j <- sort(unique(x2))

  # Vectorized binomial probability computation
  dbinom_i_mat <- outer(uniq_i, theta_seq, function(x, p) dbinom(x, N1, p))
  dbinom_j_mat <- outer(uniq_j, theta_seq, function(x, p) dbinom(x, N2, p))

  # Handle potential numerical issues
  dbinom_i_mat[is.na(dbinom_i_mat) | is.infinite(dbinom_i_mat)] <- 0
  dbinom_j_mat[is.na(dbinom_j_mat) | is.infinite(dbinom_j_mat)] <- 0

  # Map x1 and x2 to indices for efficient lookup
  i_idx <- match(x1, uniq_i)
  j_idx <- match(x2, uniq_j)

  # Compute Clopper-Pearson intervals for all points efficiently
  theta_intervals <- matrix(0, nrow = n_points, ncol = 2)
  for (k in 1:n_points) {
    x_combined <- x1[k] + x2[k]
    n_combined <- N1 + N2

    if (n_combined > 0) {
      bounds <- ClopperPearsonCI(x_combined, n_combined, gamma)
      theta_intervals[k, ] <- bounds
    } else {
      theta_intervals[k, ] <- c(0, 1)
    }
  }

  # Vectorized computation of p-values with Berger-Boos restriction
  p_values <- numeric(n_points)

  # Create a matrix to store all joint probabilities
  P_H0_all <- dbinom_i_mat[i_idx, ] * dbinom_j_mat[j_idx, ]

  # For each point, apply Berger-Boos restriction
  for (k in 1:n_points) {
    # Get Clopper-Pearson interval for this point
    theta_lower <- theta_intervals[k, 1]
    theta_upper <- theta_intervals[k, 2]

    # Find theta values within the interval
    theta_in_interval <- (theta_seq >= theta_lower) & (theta_seq <= theta_upper)

    if (any(theta_in_interval)) {
      # Extract probabilities for relevant theta values
      # For Z-pool: 1:k gives points with Z-stat >= current Z-stat (due to descending sort)
      # For Boschloo: 1:k gives points with Fisher p-val <= current Fisher p-val (due to ascending sort)
      P_H0_restricted <- P_H0_all[1:k, theta_in_interval, drop = FALSE]

      # Compute cumulative p-value efficiently
      if (ncol(P_H0_restricted) > 0) {
        # Get the probability for the current point
        current_prob <- P_H0_all[k, theta_in_interval]

        # Compute p-value as max over restricted theta range
        if (length(current_prob) > 0) {
          # For each theta in interval, compute cumulative probability
          cum_probs <- numeric(length(current_prob))
          for (t_idx in 1:length(current_prob)) {
            # Sum all probabilities >= current probability at this theta
            cum_probs[t_idx] <- sum(P_H0_restricted[, t_idx][P_H0_restricted[, t_idx] >= current_prob[t_idx]])
          }
          p_val_restricted <- max(cum_probs)
        } else {
          p_val_restricted <- 0
        }

        # Apply Berger-Boos adjustment: add gamma to the restricted p-value
        p_values[k] <- p_val_restricted + gamma
      } else {
        p_values[k] <- gamma  # If no valid theta values in interval
      }
    } else {
      # If no theta values in interval (rare case)
      p_values[k] <- gamma
    }
  }

  # Ensure p-values are in valid range [0, 1]
  p_values <- pmax(0, pmin(1, p_values))

  return(p_values)
}
