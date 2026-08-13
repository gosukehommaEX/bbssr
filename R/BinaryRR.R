#' Rejection Region for Two-Arm Trials with Binary Endpoints
#'
#' Provides a rejection region (RR) for two-arm trials with binary endpoints using various
#' exact statistical tests. Five tests are supported, each of which can be applied with a
#' one-sided or a two-sided alternative.
#'
#' @param N1 Sample size for group 1
#' @param N2 Sample size for group 2
#' @param alpha Level of significance for the alternative specified by \code{alternative}
#' @param Test Type of statistical test. Options: \code{'Chisq'}, \code{'Fisher'},
#'   \code{'Fisher-midP'}, \code{'Z-pool'}, or \code{'Boschloo'}
#' @param alternative Direction of the alternative hypothesis. Options: \code{'greater'}
#'   (default) for the one-sided alternative that the response probability of group 1
#'   exceeds that of group 2, or \code{'two.sided'}
#' @param tsmethod Convention used to construct the two-sided version of the conditional
#'   tests. Options: \code{'minlike'} (default) or \code{'central'}. Ignored when
#'   \code{alternative} is \code{'greater'}, and ignored by \code{'Chisq'} and
#'   \code{'Z-pool'}, whose two-sided versions are based on the absolute value of the Z
#'   statistic
#' @param n.grid Number of grid points used to search over the nuisance parameter of the
#'   unconditional tests. Default is 100. Ignored by the conditional tests
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure for the
#'   unconditional tests. The default of 0 disables the procedure. A positive value
#'   restricts the search over the nuisance parameter to an exact
#'   \code{100 (1 - bb.gamma)} percent confidence interval and adds \code{bb.gamma} to the
#'   resulting p-value. A common choice is 0.0001
#'
#' @return An object of class \code{bbssr_rr}, which is a logical matrix of dimension
#'   \code{(N1 + 1)} by \code{(N2 + 1)} whose entry \code{[i + 1, j + 1]} is \code{TRUE}
#'   when the null hypothesis is rejected at \code{i} responders in group 1 and \code{j}
#'   responders in group 2. The design settings are stored as attributes
#'
#' @details
#' The function supports the following five tests:
#' \itemize{
#'   \item The Pearson chi-squared test (Chisq)
#'   \item The Fisher exact test (Fisher)
#'   \item The Fisher mid-p test (Fisher-midP)
#'   \item The Z-pooled exact unconditional test (Z-pool)
#'   \item The Boschloo exact unconditional test (Boschloo)
#' }
#'
#' For the two-sided versions of the conditional tests, \code{'minlike'} sums the null
#' probabilities of all tables that are no more likely than the observed table, which is the
#' convention of \code{stats::fisher.test}, whereas \code{'central'} doubles the smaller of
#' the two one-sided tail probabilities. The two-sided versions of \code{'Chisq'} and
#' \code{'Z-pool'} order the outcomes by the absolute value of the Z statistic.
#'
#' The unconditional tests maximize the null tail probability of an ordering statistic over
#' the common response probability, which is a nuisance parameter. Outcomes sharing the same
#' value of the ordering statistic receive the same p-value.
#'
#' @examples
#' # Simple example with small sample sizes
#' RR <- BinaryRR(N1 = 5, N2 = 5, alpha = 0.025, Test = 'Chisq')
#' print(RR)
#'
#' \donttest{
#' # Two-sided Boschloo test with the Berger-Boos procedure
#' RR <- BinaryRR(N1 = 20, N2 = 10, alpha = 0.05, Test = 'Boschloo',
#'                alternative = 'two.sided', bb.gamma = 0.0001)
#' print(RR)
#' plot(RR)
#' }
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @export
#' @import fpCompare
#' @importFrom stats pnorm
BinaryRR <- function(N1, N2, alpha, Test,
                     alternative = c('greater', 'two.sided'),
                     tsmethod = c('minlike', 'central'),
                     n.grid = 100, bb.gamma = 0) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  Test <- match.arg(Test, c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo'))
  if (length(N1) != 1 || length(N2) != 1) stop('N1 and N2 must each be a single value')
  if (is.na(N1) || is.na(N2) || N1 != round(N1) || N2 != round(N2) || N1 < 1 || N2 < 1) {
    stop('N1 and N2 must be positive integers')
  }
  if (length(alpha) != 1 || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop('alpha must be a single value in (0, 1)')
  }
  if (length(n.grid) != 1 || is.na(n.grid) || n.grid < 2) {
    stop('n.grid must be a single integer of at least 2')
  }
  if (length(bb.gamma) != 1 || is.na(bb.gamma) || bb.gamma < 0 || bb.gamma >= alpha) {
    stop('bb.gamma must be a single value satisfying 0 <= bb.gamma < alpha')
  }
  unconditional <- Test %in% c('Z-pool', 'Boschloo')
  if (bb.gamma > 0 && !unconditional) {
    warning('bb.gamma applies only to the Z-pool and Boschloo tests and is ignored here')
  }
  N1 <- as.integer(N1)
  N2 <- as.integer(N2)
  n.grid <- as.integer(n.grid)
  if (Test == 'Chisq') {
    Z <- zstat(N1, N2)
    p.val <- if (alternative == 'greater') {
      pnorm(Z, lower.tail = FALSE)
    } else {
      pmin(2 * pnorm(abs(Z), lower.tail = FALSE), 1)
    }
  } else if (Test == 'Fisher') {
    p.val <- fisher_pvalue(N1, N2, alternative, tsmethod, midp = FALSE)
  } else if (Test == 'Fisher-midP') {
    p.val <- fisher_pvalue(N1, N2, alternative, tsmethod, midp = TRUE)
  } else if (Test == 'Z-pool') {
    Z <- zstat(N1, N2)
    stat <- if (alternative == 'greater') Z else abs(Z)
    p.val <- unconditional_pvalue(stat, N1, N2, n.grid, bb.gamma, decreasing = TRUE)
  } else {
    stat <- fisher_pvalue(N1, N2, alternative, tsmethod, midp = FALSE)
    p.val <- unconditional_pvalue(stat, N1, N2, n.grid, bb.gamma, decreasing = FALSE)
  }
  RR <- (p.val %<<% alpha)
  dim(RR) <- c(N1 + 1L, N2 + 1L)
  dimnames(RR) <- list(x1 = as.character(0:N1), x2 = as.character(0:N2))
  attr(RR, 'N1') <- N1
  attr(RR, 'N2') <- N2
  attr(RR, 'alpha') <- alpha
  attr(RR, 'Test') <- Test
  attr(RR, 'alternative') <- alternative
  attr(RR, 'tsmethod') <- tsmethod
  attr(RR, 'n.grid') <- n.grid
  attr(RR, 'bb.gamma') <- bb.gamma
  attr(RR, 'p.value') <- p.val
  class(RR) <- c('bbssr_rr', 'matrix', 'array')
  RR
}
