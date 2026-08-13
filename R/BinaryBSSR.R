#' Sample Size Re-estimation from Observed Blinded Interim Data
#'
#' Re-estimates the sample size of an ongoing two-arm trial with a binary endpoint from the
#' blinded data available at an interim analysis, and reports how many patients still have
#' to be enrolled in each group during the second stage. Only the total number of patients
#' and the total number of responders are required, so the treatment allocation remains
#' concealed.
#'
#' @param n1 Number of patients of group 1 observed at the interim analysis
#' @param n2 Number of patients of group 2 observed at the interim analysis
#' @param S Total number of responders observed at the interim analysis, pooled over both
#'   groups
#' @param Delta.A Assumed treatment effect (risk difference) used to split the blinded
#'   pooled proportion into group-specific proportions
#' @param r Allocation ratio to group 1 (i.e., allocation ratio of group 1:group 2 = r:1,
#'   r > 0)
#' @param alpha Level of significance for the alternative specified by \code{alternative}
#' @param tar.power Target power
#' @param Test Type of statistical test. Options: \code{'Chisq'}, \code{'Fisher'},
#'   \code{'Fisher-midP'}, \code{'Z-pool'}, or \code{'Boschloo'}
#' @param restricted Logical. If \code{TRUE}, the re-estimated sample size is not allowed
#'   to fall below the planned sample size given by \code{N1} and \code{N2}. Default is
#'   \code{FALSE}
#' @param N1 Planned sample size of group 1. Required when \code{restricted} is \code{TRUE}
#' @param N2 Planned sample size of group 2. Required when \code{restricted} is \code{TRUE}
#' @param alternative Direction of the alternative hypothesis. Options: \code{'greater'}
#'   (default) or \code{'two.sided'}
#' @param tsmethod Convention used to construct the two-sided version of the conditional
#'   tests. Options: \code{'minlike'} (default) or \code{'central'}
#' @param n.grid Number of grid points used to search over the nuisance parameter of the
#'   unconditional tests. Default is 100
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure. The default of
#'   0 disables the procedure
#'
#' @return An object of class \code{bbssr_bssr}, a data frame with one row containing:
#' \describe{
#'   \item{n1}{Interim sample size of group 1}
#'   \item{n2}{Interim sample size of group 2}
#'   \item{n}{Total interim sample size}
#'   \item{S}{Total number of interim responders}
#'   \item{hat.p}{Blinded estimate of the pooled response probability}
#'   \item{hat.p1}{Recovered response probability of group 1}
#'   \item{hat.p2}{Recovered response probability of group 2}
#'   \item{N1.re}{Re-estimated total sample size of group 1}
#'   \item{N2.re}{Re-estimated total sample size of group 2}
#'   \item{N.re}{Re-estimated total sample size}
#'   \item{n1.stage2}{Number of additional patients to enrol in group 1}
#'   \item{n2.stage2}{Number of additional patients to enrol in group 2}
#'   \item{n.stage2}{Total number of additional patients to enrol}
#'   \item{N1.final}{Final sample size of group 1}
#'   \item{N2.final}{Final sample size of group 2}
#'   \item{N.final}{Final total sample size}
#'   \item{Power}{Exact power at the final sample size under the recovered proportions}
#' }
#'
#' @details
#' The blinded estimate of the pooled response probability is \code{hat.p = S / (n1 + n2)}.
#' Group-specific proportions are recovered as
#' \code{hat.p1 = hat.p + Delta.A / (1 + r)} and
#' \code{hat.p2 = hat.p - r Delta.A / (1 + r)}, truncated to the unit interval, and the
#' sample size is then re-estimated by \code{\link{BinarySampleSize}}.
#'
#' Under the unrestricted rule the second-stage sample size is the re-estimated sample size
#' minus what has already been observed, and is never negative. Under the restricted rule
#' the re-estimated sample size is first raised to the planned sample size, so the trial can
#' only grow.
#'
#' While \code{\link{BinaryPowerBSSR}} evaluates the operating characteristics of a BSSR
#' design at the planning stage, this function is applied once, to the data of a trial that
#' is under way.
#'
#' @examples
#' # Interim data: 20 patients per group, 11 responders in total
#' BinaryBSSR(n1 = 20, n2 = 20, S = 11, Delta.A = 0.3, r = 1,
#'            alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
#'
#' \donttest{
#' # Restricted rule with a planned sample size of 40 per group
#' BinaryBSSR(n1 = 20, n2 = 20, S = 11, Delta.A = 0.3, r = 1,
#'            alpha = 0.025, tar.power = 0.8, Test = 'Boschloo',
#'            restricted = TRUE, N1 = 40, N2 = 40)
#' }
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @seealso \code{\link{BinaryPowerBSSR}}, \code{\link{BinarySampleSize}}
#' @export
BinaryBSSR <- function(n1, n2, S, Delta.A, r, alpha, tar.power, Test,
                       restricted = FALSE, N1 = NULL, N2 = NULL,
                       alternative = c('greater', 'two.sided'),
                       tsmethod = c('minlike', 'central'),
                       n.grid = 100, bb.gamma = 0) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  if (length(n1) != 1 || length(n2) != 1 || length(S) != 1) {
    stop('n1, n2 and S must each be a single value')
  }
  if (n1 != round(n1) || n2 != round(n2) || n1 < 1 || n2 < 1) {
    stop('n1 and n2 must be positive integers')
  }
  if (S != round(S) || S < 0 || S > n1 + n2) {
    stop('S must be an integer between 0 and n1 + n2')
  }
  if (length(Delta.A) != 1 || Delta.A <= 0 || Delta.A >= 1) {
    stop('Delta.A must be a single value in (0, 1)')
  }
  if (restricted && (is.null(N1) || is.null(N2))) {
    stop('N1 and N2 must be supplied when restricted is TRUE')
  }
  n1 <- as.integer(n1)
  n2 <- as.integer(n2)
  S <- as.integer(S)
  # Blinded estimate of the pooled response probability
  hat.p <- S / (n1 + n2)
  hat.p1 <- pmin(1, hat.p + (1 / (1 + r)) * Delta.A)
  hat.p2 <- pmax(0, hat.p - (r / (1 + r)) * Delta.A)
  # Sample size re-estimation
  ss <- BinarySampleSize(hat.p1, hat.p2, r, alpha, tar.power, Test,
                         alternative, tsmethod, n.grid, bb.gamma)
  N1.re <- ss[['N1']]
  N2.re <- ss[['N2']]
  # Second-stage sample sizes
  n2.stage2 <- if (restricted) as.integer(max(N2, N2.re) - n2) else as.integer(max(n2, N2.re) - n2)
  n1.stage2 <- as.integer(ceiling(r * n2.stage2))
  N1.final <- n1 + n1.stage2
  N2.final <- n2 + n2.stage2
  Power <- BinaryPower(hat.p1, hat.p2, N1.final, N2.final, alpha, Test,
                       alternative, tsmethod, n.grid, bb.gamma)$Power
  out <- data.frame(
    n1 = n1, n2 = n2, n = n1 + n2, S = S,
    hat.p = hat.p, hat.p1 = hat.p1, hat.p2 = hat.p2,
    N1.re = N1.re, N2.re = N2.re, N.re = N1.re + N2.re,
    n1.stage2 = n1.stage2, n2.stage2 = n2.stage2, n.stage2 = n1.stage2 + n2.stage2,
    N1.final = N1.final, N2.final = N2.final, N.final = N1.final + N2.final,
    Power = Power,
    stringsAsFactors = FALSE
  )
  attr(out, 'Test') <- Test
  attr(out, 'alternative') <- alternative
  attr(out, 'alpha') <- alpha
  attr(out, 'tar.power') <- tar.power
  attr(out, 'restricted') <- restricted
  attr(out, 'Delta.A') <- Delta.A
  attr(out, 'r') <- r
  class(out) <- c('bbssr_bssr', 'data.frame')
  out
}
