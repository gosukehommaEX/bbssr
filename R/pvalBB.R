#' P-value Computation with Berger-Boos Approach
#'
#' Computes p-values using the Berger-Boos approach for exact unconditional tests
#' with advanced optimization techniques for maximum performance. This function
#' provides significant computational efficiency while maintaining statistical
#' accuracy for both Z-pooled and Boschloo unconditional exact tests.
#'
#' @param x1 Numeric vector of success counts for group 1 (pre-sorted by test statistic)
#' @param x2 Numeric vector of success counts for group 2 (pre-sorted by test statistic)
#' @param N1 Integer, sample size for group 1
#' @param N2 Integer, sample size for group 2
#' @param gamma Numeric, Berger-Boos parameter for confidence interval adjustment (0 < gamma < 1)
#'
#' @return Numeric vector of p-values corresponding to each (x1, x2) pair
#'
#' @details
#' This function implements the Berger-Boos approach for computing exact unconditional
#' p-values with the following key optimizations:
#' \itemize{
#'   \item Vectorized Clopper-Pearson confidence interval computation
#'   \item Pre-computed probability lookup tables for all parameter combinations
#'   \item Matrix-based theta parameter filtering for efficient restriction
#'   \item Advanced memory access pattern optimization
#'   \item Elimination of redundant binomial probability calculations
#' }
#'
#' The input vectors x1 and x2 must be pre-sorted according to the test statistic:
#' \itemize{
#'   \item For Z-pooled test: sorted by Z-statistic in descending order
#'   \item For Boschloo test: sorted by Fisher p-values in ascending order
#' }
#'
#' The function uses Clopper-Pearson confidence intervals to restrict the nuisance
#' parameter space, then computes the maximum p-value over this restricted range.
#' The final p-value includes the gamma adjustment per Berger-Boos methodology.
#'
#' @examples
#' \dontrun{
#' # Example for Z-pooled test (sorted by Z-statistic descending)
#' x1.vals <- c(12, 10, 8, 6)
#' x2.vals <- c(5, 4, 3, 2)
#' p.vals <- pvalBB(x1.vals, x2.vals, N1 = 25, N2 = 20, gamma = 0.0001)
#'
#' # Example for Boschloo test (sorted by Fisher p-value ascending)
#' x1.vals <- c(6, 8, 10, 12)
#' x2.vals <- c(2, 3, 4, 5)
#' p.vals <- pvalBB(x1.vals, x2.vals, N1 = 25, N2 = 20, gamma = 0.0001)
#' }
#'
#' @references
#' Berger, R. L., & Boos, D. D. (1994). P values maximized over a confidence
#' set for the nuisance parameter. Journal of the American Statistical
#' Association, 89(427), 1012-1016.
#'
#' Clopper, C., & Pearson, E. S. (1934). The use of confidence or fiducial
#' limits illustrated in the case of the binomial. Biometrika, 26(4), 404-413.
#'
#' @export
#' @importFrom stats qbeta dbinom
pvalBB <- function(x1, x2, N1, N2, gamma) {

  # Input validation with informative error messages
  if (length(x1) != length(x2)) {
    stop("x1 and x2 must have the same length", call. = FALSE)
  }

  n.pts <- length(x1)
  if (n.pts == 0) return(numeric(0))

  if (any(x1 < 0, x2 < 0, na.rm = TRUE)) {
    stop("All values in x1 and x2 must be non-negative", call. = FALSE)
  }
  if (any(x1 > N1, x2 > N2, na.rm = TRUE)) {
    stop("Values in x1 and x2 cannot exceed their respective sample sizes", call. = FALSE)
  }
  if (gamma <= 0 || gamma >= 1) {
    stop("gamma must be between 0 and 1", call. = FALSE)
  }

  # Create theta parameter grid
  theta.seq <- seq(0, 1, length.out = 100)

  # Vectorized Clopper-Pearson confidence interval computation
  x.comb <- x1 + x2
  n.comb <- N1 + N2

  # Compute lower bounds efficiently
  theta.lower <- numeric(n.pts)
  non.zero.mask <- x.comb > 0
  if (any(non.zero.mask)) {
    theta.lower[non.zero.mask] <- qbeta(
      gamma/2, x.comb[non.zero.mask],
      n.comb - x.comb[non.zero.mask] + 1
    )
  }
  # theta.lower[!non.zero.mask] remains 0 by default

  # Compute upper bounds efficiently
  theta.upper <- rep(1, n.pts)
  non.full.mask <- x.comb < n.comb
  if (any(non.full.mask)) {
    theta.upper[non.full.mask] <- qbeta(
      1 - gamma/2, x.comb[non.full.mask] + 1,
      n.comb - x.comb[non.full.mask]
    )
  }

  # Vectorized computation of all joint probabilities
  # Create expanded parameter vectors for single dbinom call
  x1.exp <- rep(x1, each = 100)
  x2.exp <- rep(x2, each = 100)
  theta.exp <- rep(theta.seq, times = n.pts)

  # Compute all binomial probabilities at once
  prob1.vec <- dbinom(x1.exp, N1, theta.exp)
  prob2.vec <- dbinom(x2.exp, N2, theta.exp)
  joint.prob.vec <- prob1.vec * prob2.vec

  # Handle numerical issues (inf, -inf, NaN)
  joint.prob.vec[!is.finite(joint.prob.vec)] <- 0

  # Reshape into matrix for efficient indexing
  joint.probs <- matrix(joint.prob.vec, nrow = n.pts, ncol = 100, byrow = TRUE)

  # Create theta validity mask using matrix operations
  theta.mat <- matrix(theta.seq, nrow = n.pts, ncol = 100, byrow = TRUE)
  lower.mat <- matrix(theta.lower, nrow = n.pts, ncol = 100)
  upper.mat <- matrix(theta.upper, nrow = n.pts, ncol = 100)

  valid.theta.mask <- (theta.mat >= lower.mat) & (theta.mat <= upper.mat)

  # Optimized p-value computation
  p.vals <- numeric(n.pts)

  for (k in 1:n.pts) {
    # Find valid theta indices for current point
    valid.cols <- which(valid.theta.mask[k, ])

    if (length(valid.cols) == 0) {
      p.vals[k] <- gamma
      next
    }

    # Compute maximum p-value over valid theta range
    max.pval <- 0

    for (col in valid.cols) {
      curr.prob <- joint.probs[k, col]

      # Find all probabilities >= current probability in relevant range
      relevant.probs <- joint.probs[1:k, col]
      cum.prob <- sum(relevant.probs[relevant.probs >= curr.prob])

      max.pval <- max(max.pval, cum.prob)
    }

    # Apply Berger-Boos adjustment
    p.vals[k] <- max.pval + gamma
  }

  # Ensure p-values are within valid range [0, 1]
  p.vals <- pmax(0, pmin(1, p.vals))

  return(p.vals)
}
