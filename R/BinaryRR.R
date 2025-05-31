#' Rejection Region for Two-Arm Trials with Binary Endpoints
#'
#' Provides a rejection region (RR) for two-arm trials with binary endpoints using
#' various exact statistical tests. The function supports five different one-sided tests
#' with optional Berger-Boos approach based on Clopper-Pearson confidence intervals.
#'
#' @param N1 Sample size for group 1.
#' @param N2 Sample size for group 2.
#' @param alpha One-sided level of significance.
#' @param Test Type of statistical test. Options: 'Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', or 'Boschloo'.
#' @param BB Logical. If TRUE, apply Berger-Boos approach for 'Z-pool' and 'Boschloo' tests (default: FALSE).
#' @param gamma Berger-Boos parameter for confidence interval adjustment (default: 0.0001).
#'
#' @return A logical matrix representing the rejection region (RR). Matrix dimensions
#' are (N1+1) x (N2+1), where TRUE indicates rejection of the null hypothesis.
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
#' When BB = TRUE, the Berger-Boos approach based on Clopper-Pearson
#' confidence intervals is applied to 'Z-pool' and 'Boschloo' tests for improved
#' computational efficiency while maintaining statistical validity.
#'
#' @examples
#' \dontrun{
#' # Calculate rejection region for Boschloo test
#' N1 <- 20
#' N2 <- 10
#' alpha <- 0.025
#' Test <- 'Boschloo'
#' RR <- BinaryRR(N1, N2, alpha, Test)
#' print(RR)
#'
#' # Calculate rejection region with Berger-Boos approach
#' RR_BB <- BinaryRR(N1, N2, alpha, Test, BB = TRUE)
#' print(RR_BB)
#' }
#'
#' @export
#' @import fpCompare
#' @import stats pnorm dbinom phyper dhyper qbeta
BinaryRR <- function(N1, N2, alpha, Test, BB = FALSE, gamma = 0.0001) {
  # Input validation
  if (N1 <= 0 || N2 <= 0) stop('Sample sizes N1 and N2 must be positive')
  if (alpha <= 0 || alpha >= 1) stop('alpha must be between 0 and 1')
  if (!Test %in% c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')) {
    stop("Test must be one of: 'Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', or 'Boschloo'")
  }
  if (gamma <= 0 || gamma >= 1) stop('gamma must be between 0 and 1')

  if((Test == 'Chisq') | (Test == 'Z-pool')) {
    # Test statistics for the Chisq over all combinations of i(=0,...,N1) and j(=0,...,N2)
    Z.ij <- '/'(
      outer(N2 * (0:N1), N1 * (0:N2), '-') / (N1 * N2),
      sqrt(outer(0:N1, 0:N2, '+') / (N1 * N2) * (1 - outer(0:N1, 0:N2, '+') / (N1 + N2)))
    )
    Z.ij[is.na(Z.ij)] <- 0

    if(Test == 'Chisq') {
      # p-values for the Chisq
      p.val <- pnorm(Z.ij, lower.tail = FALSE)
    } else if(Test == 'Z-pool') {
      # Since zero and negative values of test statistics must not be statistically significant, they are omitted
      Z.ij.posi <- Z.ij[Z.ij %>>% 0]
      order.Z.ij.posi <- order(Z.ij.posi, decreasing = TRUE)
      i <- (row(Z.ij)[Z.ij %>>% 0] - 1)[order.Z.ij.posi]
      j <- (col(Z.ij)[Z.ij %>>% 0] - 1)[order.Z.ij.posi]

      if (BB) {
        # Apply Berger-Boos approach for Z-pool test
        p.ij <- pvalBB(i, j, N1, N2, gamma)
      } else {
        # Original Z-pool computation
        # Calculate P_H0(X1 = i, X2 = j | theta)
        uniq.i <- sort(unique(i))
        dbinom.i <- sapply(seq(0, 1, l = 100), function(theta) dbinom(uniq.i, N1, theta))
        uniq.j <- sort(unique(j))
        dbinom.j <- sapply(seq(0, 1, l = 100), function(theta) dbinom(uniq.j, N2, theta))
        P_H0 <- dbinom.i[i - min(uniq.i) + 1, ] * dbinom.j[j + 1, ]
        # p-values for all possible values of theta
        p.ij <- apply(apply(P_H0, 2, cumsum), 1, max)
      }

      # p-values for the Z-pool
      p.val <- 1 ^ Z.ij
      p.val[cbind(i + 1, j + 1)] <- p.ij
    }
  } else if(Test == 'Fisher-midP') {
    # p-values for the Fisher-midP
    p.val <- outer(0:N1, 0:N2, function(i, j) phyper(i, N1, N2, i + j, lower.tail = FALSE) + 0.5 * dhyper(i, N1, N2, i + j))
  } else if((Test == 'Fisher') | (Test == 'Boschloo')) {
    # p-values for the Fisher over all combinations of i(=0,...,N1) and j(=0,...,N2)
    p.fisher <- outer(0:N1, 0:N2, function(i, j) phyper(i - 1, N1, N2, i + j, lower.tail = FALSE))

    if(Test == 'Fisher') {
      # p-values for the Fisher
      p.val <- p.fisher
    } else if(Test == 'Boschloo') {
      # Since p-values satisfying hat{p}_{1} - hat{p}_{2} <= 0 must not be statistically significant, they are omitted
      p.max.boschloo <- min(p.fisher[(outer(N2 * (0:N1), N1 * (0:N2), '-') / (N1 * N2)) %<=% 0])
      p.fisher.posi <- p.fisher[p.fisher %<<% p.max.boschloo]
      order.p.fisher.posi <- order(p.fisher.posi, decreasing = FALSE)
      i <- (row(p.fisher)[p.fisher %<<% p.max.boschloo] - 1)[order.p.fisher.posi]
      j <- (col(p.fisher)[p.fisher %<<% p.max.boschloo] - 1)[order.p.fisher.posi]

      if (BB) {
        # Apply Berger-Boos approach for Boschloo test
        p.ij <- pvalBB(i, j, N1, N2, gamma)
      } else {
        # Original Boschloo computation
        # Calculate P_H0(X1 = i, X2 = j | theta)
        uniq.i <- sort(unique(i))
        dbinom.i <- sapply(seq(0, 1, l = 100), function(theta) dbinom(uniq.i, N1, theta))
        uniq.j <- sort(unique(j))
        dbinom.j <- sapply(seq(0, 1, l = 100), function(theta) dbinom(uniq.j, N2, theta))
        P_H0 <- dbinom.i[i - min(uniq.i) + 1, ] * dbinom.j[j + 1, ]
        # p-values for all possible values of theta
        p.ij <- apply(apply(P_H0, 2, cumsum), 1, max)
      }

      # p-values for the Boschloo
      p.val <- 1 ^ p.fisher
      p.val[cbind(i + 1, j + 1)] <- p.ij
    }
  }

  # Rejection region
  RR <- (p.val %<<% alpha)

  # Return RR
  return(RR)
}
