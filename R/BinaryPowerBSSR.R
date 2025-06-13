#' Power Calculation for Two-Arm Trials with Binary Endpoints Using Blinded Sample Size Re-estimation (BSSR)
#'
#' Calculates the power for two-arm trials with binary endpoints when blinded sample size
#' re-estimation (BSSR) is implemented. The function supports five different statistical tests
#' and allows for both restricted and unrestricted designs with optional weighted approaches.
#'
#' @param asmd.p1 Assumed proportion of responders for group 1
#' @param asmd.p2 Assumed proportion of responders for group 2
#' @param p Vector of pooled proportions of responders from both groups (can specify multiple values)
#' @param Delta.A Assumed treatment effect (risk difference)
#' @param Delta.T True treatment effect (risk difference)
#' @param N1 Initial sample size of group 1
#' @param N2 Initial sample size of group 2
#' @param omega Fraction of sample size used for interim analysis (i.e., for BSSR)
#' @param r Allocation ratio to group 1
#' @param alpha One-sided level of significance
#' @param tar.power Target power
#' @param Test Type of statistical test. Options: 'Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', or 'Boschloo'
#' @param restricted Logical. If TRUE, restricted design is chosen
#' @param weighted Logical. If TRUE, weighted approach is chosen
#'
#' @return A data frame containing:
#' \describe{
#'   \item{p1}{True probability of responders for group 1}
#'   \item{p2}{True probability of responders for group 2}
#'   \item{p}{True probability of pooled responders from both groups}
#'   \item{power.BSSR}{Power for BSSR design}
#'   \item{power.TRAD}{Power for traditional design}
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
#' @examples
#' \dontrun{
#' result <- BinaryPowerBSSR(
#'   asmd.p1 = 0.45,
#'   asmd.p2 = 0.09,
#'   p = seq(0.14, 0.23, by = 0.01),
#'   Delta.A = 0.36,
#'   Delta.T = 0.36,
#'   N1 = 24,
#'   N2 = 24,
#'   omega = 0.5,
#'   r = 1,
#'   alpha = 0.025,
#'   tar.power = 0.8,
#'   Test = 'Z-pool',
#'   restricted = FALSE,
#'   weighted = TRUE
#' )
#' print(result)
#' }
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @export
#' @import fpCompare
#' @import stats dbinom pbinom
BinaryPowerBSSR <- function(asmd.p1, asmd.p2, p, Delta.A, Delta.T, N1, N2, omega, r, alpha, tar.power, Test, restricted, weighted) {
  # Initial sample sizes used for BSSR
  N11 <- ceiling(omega * N1)
  N12 <- ceiling(omega * N2)
  # Estimate pooled proportion under blinded fashion
  hat.p <- c(outer(0:N11, 0:N12, '+') / (N11 + N12))
  # Derive hat{p}_{1} and hat{p}_{2} using hat{p} and Delta
  hat.p1 <- pmin(1, hat.p + (1 / (1 + r)) * Delta.A)
  hat.p2 <- pmax(0, hat.p - (r / (1 + r)) * Delta.A)
  # Sample size re-estimation
  BSSR.N2 <- sapply(seq(hat.p), function(i) BinarySampleSize(hat.p1[i], hat.p2[i], r, alpha, tar.power, Test)[['N2']])
  if(restricted == TRUE) { N22 <- pmax(N2, BSSR.N2) - N12 } else { N22 <- pmax(N12, BSSR.N2) - N12 }
  N21 <- r * N22
  # Final sample sizes
  hat.N1 <- N11 + N21
  hat.N2 <- N12 + N22
  hat.N <- hat.N1 + hat.N2
  if(weighted == TRUE) {
    # Define weights
    w.h <- c(outer(dbinom(0:N11, N11, asmd.p1), dbinom(0:N12, N12, asmd.p2)))
    N.w <- ceiling(sum(hat.N * w.h))
    # Final sample sizes
    hat.N <- pmax(hat.N, N.w)
    hat.N2 <- ceiling(hat.N / (1 + r))
    hat.N1 <- ceiling(r * hat.N2)
    N21 <- hat.N1 - N11
    N22 <- hat.N2 - N12
  }
  # True proportion for each treatment group
  p1 <- p + (1 / (1 + r)) * Delta.T
  p2 <- p - (r / (1 + r)) * Delta.T
  # Omit scenarios with negative or greater than 1 p_{j} values
  omit.ID <- union(which(p1 %>>% 1), which(p2 %<<% 0))
  if(length(omit.ID) %!=% 0) {
    p1 <- p1[-omit.ID]
    p2 <- p2[-omit.ID]
    p <- p[-omit.ID]
  }
  # Set probability mass functions of binomial distributions for the interim analysis
  dbinom1 <- outer(X = 0:N11, Y = p1, function(X, Y) dbinom(X, N11, Y))
  dbinom2 <- outer(X = 0:N12, Y = p2, function(X, Y) dbinom(X, N12, Y))
  # Powers given \hat{N}_{1} and \hat{N}_{2}
  x11 <- c(row(matrix(0, nrow = N11 + 1, ncol = N12 + 1)) - 1)
  x12 <- c(col(matrix(0, nrow = N11 + 1, ncol = N12 + 1)) - 1)
  power.stage2 <- do.call(
    cbind,
    lapply(seq(hat.p), function(i) {
      # Set rejection region
      RR <- BinaryRR(hat.N1[i], hat.N2[i], alpha, Test)
      RR <- as.matrix(RR[(x11[i] + 1):(x11[i] + 1 + N21[i]), (x12[i] + 1):(x12[i] + 1 + N22[i])])
      # Return power
      sapply(seq(p), function(j) {
        c(dbinom(0:N21[i], N21[i], p1[j]) %*% pbinom(rowSums(RR) - 1, N22[i], p2[j]))
      })
    })
  )
  # Final power for BSSR
  power.BSSR <- sapply(seq(p), function(k) {
    sum(c(dbinom1[, k] %o% dbinom2[, k]) * power.stage2[k, ], na.rm = TRUE)
  })
  # Powers without BSSR (i.e., traditional design)
  power.TRAD <- BinaryPower(p1, p2, N1, N2, alpha, Test)
  # Output
  out <- data.frame(p1, p2, p, power.BSSR, power.TRAD)
  return(out)
}
