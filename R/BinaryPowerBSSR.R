#' Power of a Blinded Sample Size Re-estimation Design for Binary Endpoints
#'
#' Calculates the power of a two-arm trial with a binary endpoint when blinded sample size
#' re-estimation (BSSR) is implemented, together with the power of the corresponding
#' fixed-sample design. Five exact tests are supported, each of which can be applied with a
#' one-sided or a two-sided alternative, under either a restricted or an unrestricted
#' design rule.
#'
#' @param asmd.p1 Assumed proportion of responders for group 1
#' @param asmd.p2 Assumed proportion of responders for group 2
#' @param p Vector of true pooled proportions of responders from both groups
#' @param Delta.A Assumed treatment effect (risk difference)
#' @param Delta.T True treatment effect (risk difference)
#' @param N1 Initial sample size of group 1
#' @param N2 Initial sample size of group 2
#' @param omega Fraction of the initial sample size observed at the interim analysis
#' @param r Allocation ratio to group 1
#' @param alpha Level of significance for the alternative specified by \code{alternative}
#' @param tar.power Target power
#' @param Test Type of statistical test. Options: \code{'Chisq'}, \code{'Fisher'},
#'   \code{'Fisher-midP'}, \code{'Z-pool'}, or \code{'Boschloo'}
#' @param restricted Logical. If \code{TRUE}, the re-estimated sample size is not allowed
#'   to fall below the initial sample size. Default is \code{FALSE}
#' @param alternative Direction of the alternative hypothesis. Options: \code{'greater'}
#'   (default) or \code{'two.sided'}
#' @param tsmethod Convention used to construct the two-sided version of the conditional
#'   tests. Options: \code{'minlike'} (default) or \code{'central'}
#' @param n.grid Number of grid points used to search over the nuisance parameter of the
#'   unconditional tests. Default is 100
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure. The default of
#'   0 disables the procedure
#'
#' @return An object of class \code{bbssr_powerbssr}, a data frame with one row per element
#'   of \code{p} containing:
#' \describe{
#'   \item{p1}{True probability of responders for group 1}
#'   \item{p2}{True probability of responders for group 2}
#'   \item{p}{True pooled probability of responders from both groups}
#'   \item{power.BSSR}{Power of the BSSR design}
#'   \item{power.TRAD}{Power of the fixed-sample design}
#'   \item{E.N}{Expected total sample size of the BSSR design}
#' }
#'
#' @details
#' At the interim analysis the pooled number of responders is observed without unblinding.
#' The pooled proportion is combined with the assumed treatment effect \code{Delta.A} to
#' recover group-specific proportions, from which the sample size is re-estimated. The power
#' is then averaged over the distribution of the interim outcome.
#'
#' Setting \code{Delta.T} to 0 gives the exact type I error rate of the design.
#'
#' The weighted approach available in earlier versions of the package has been removed.
#'
#' @examples
#' # Small BSSR calculation with the chi-squared test
#' BinaryPowerBSSR(
#'   asmd.p1 = 0.6, asmd.p2 = 0.3, p = 0.45,
#'   Delta.A = 0.3, Delta.T = 0.3,
#'   N1 = 5, N2 = 5, omega = 0.5, r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
#' )
#'
#' \donttest{
#' res <- BinaryPowerBSSR(
#'   asmd.p1 = 0.45, asmd.p2 = 0.09, p = seq(0.14, 0.23, by = 0.01),
#'   Delta.A = 0.36, Delta.T = 0.36,
#'   N1 = 24, N2 = 24, omega = 0.5, r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Z-pool'
#' )
#' print(res)
#' plot(res)
#' }
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @export
#' @import fpCompare
#' @importFrom stats dbinom
BinaryPowerBSSR <- function(asmd.p1, asmd.p2, p, Delta.A, Delta.T, N1, N2, omega, r,
                            alpha, tar.power, Test, restricted = FALSE,
                            alternative = c('greater', 'two.sided'),
                            tsmethod = c('minlike', 'central'),
                            n.grid = 100, bb.gamma = 0) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  if (length(omega) != 1 || omega <= 0 || omega > 1) {
    stop('omega must be a single value in (0, 1]')
  }
  # Initial sample sizes used for BSSR
  N11 <- as.integer(ceiling(omega * N1))
  N12 <- as.integer(ceiling(omega * N2))
  # Estimate the pooled proportion in a blinded fashion
  hat.p <- c(outer(0:N11, 0:N12, '+') / (N11 + N12))
  # Derive hat{p}_{1} and hat{p}_{2} using hat{p} and Delta.A
  hat.p1 <- pmin(1, hat.p + (1 / (1 + r)) * Delta.A)
  hat.p2 <- pmax(0, hat.p - (r / (1 + r)) * Delta.A)
  # Sample size re-estimation, evaluated once per distinct pair of interim proportions
  key.p <- paste(hat.p1, hat.p2, sep = '_')
  uniq.p <- !duplicated(key.p)
  ss.cache <- new.env(parent = emptyenv())
  for (u in which(uniq.p)) {
    assign(
      key.p[u],
      BinarySampleSize(hat.p1[u], hat.p2[u], r, alpha, tar.power, Test,
                       alternative, tsmethod, n.grid, bb.gamma)[['N2']],
      envir = ss.cache
    )
  }
  BSSR.N2 <- vapply(key.p, function(k) as.numeric(get(k, envir = ss.cache)),
                    numeric(1), USE.NAMES = FALSE)
  if (restricted == TRUE) N22 <- pmax(N2, BSSR.N2) - N12 else N22 <- pmax(N12, BSSR.N2) - N12
  N22 <- as.integer(N22)
  N21 <- as.integer(ceiling(r * N22))
  # Final sample sizes
  hat.N1 <- N11 + N21
  hat.N2 <- N12 + N22
  hat.N <- hat.N1 + hat.N2
  # True proportion for each treatment group
  p1 <- p + (1 / (1 + r)) * Delta.T
  p2 <- p - (r / (1 + r)) * Delta.T
  # Omit scenarios with p_{j} outside the unit interval
  omit.ID <- union(which(p1 %>>% 1), which(p2 %<<% 0))
  if (length(omit.ID) %!=% 0) {
    p1 <- p1[-omit.ID]
    p2 <- p2[-omit.ID]
    p <- p[-omit.ID]
  }
  if (length(p) == 0) stop('no scenario has both p1 and p2 inside the unit interval')
  # Probability mass functions of the interim outcome
  dbinom1 <- outer(X = 0:N11, Y = p1, function(X, Y) dbinom(X, N11, Y))
  dbinom2 <- outer(X = 0:N12, Y = p2, function(X, Y) dbinom(X, N12, Y))
  # Interim responder counts, in the same order as hat.p
  x11 <- c(row(matrix(0L, nrow = N11 + 1L, ncol = N12 + 1L)) - 1L)
  x12 <- c(col(matrix(0L, nrow = N11 + 1L, ncol = N12 + 1L)) - 1L)
  # Rejection regions, evaluated once per distinct pair of final sample sizes
  key.N <- paste(hat.N1, hat.N2, sep = '_')
  uniq.N <- !duplicated(key.N)
  rr.cache <- new.env(parent = emptyenv())
  for (u in which(uniq.N)) {
    RR <- BinaryRR(hat.N1[u], hat.N2[u], alpha, Test, alternative, tsmethod, n.grid, bb.gamma)
    assign(key.N[u], matrix(as.vector(RR), nrow = hat.N1[u] + 1L, ncol = hat.N2[u] + 1L),
           envir = rr.cache)
  }
  # Conditional power of the second stage given each interim outcome
  power.stage2 <- do.call(
    cbind,
    lapply(seq_along(hat.p), function(i) {
      rr.full <- get(key.N[i], envir = rr.cache)
      rows <- (x11[i] + 1L):(x11[i] + 1L + N21[i])
      cols <- (x12[i] + 1L):(x12[i] + 1L + N22[i])
      rr <- matrix(rr.full[rows, cols], nrow = length(rows), ncol = length(cols))
      vapply(
        seq_along(p),
        function(j) power_from_rr(rr, dbinom(0:N21[i], N21[i], p1[j]),
                                  dbinom(0:N22[i], N22[i], p2[j])),
        numeric(1)
      )
    })
  )
  # Average over the distribution of the interim outcome
  power.BSSR <- vapply(seq_along(p), function(k) {
    sum(c(dbinom1[, k] %o% dbinom2[, k]) * power.stage2[k, ], na.rm = TRUE)
  }, numeric(1))
  E.N <- vapply(seq_along(p), function(k) {
    sum(c(dbinom1[, k] %o% dbinom2[, k]) * hat.N)
  }, numeric(1))
  # Power of the fixed-sample design
  power.TRAD <- BinaryPower(p1, p2, N1, N2, alpha, Test, alternative, tsmethod,
                            n.grid, bb.gamma)$Power
  out <- data.frame(p1, p2, p, power.BSSR, power.TRAD, E.N)
  attr(out, 'Test') <- Test
  attr(out, 'alternative') <- alternative
  attr(out, 'alpha') <- alpha
  attr(out, 'tar.power') <- tar.power
  attr(out, 'restricted') <- restricted
  attr(out, 'N1') <- N1
  attr(out, 'N2') <- N2
  attr(out, 'omega') <- omega
  attr(out, 'Delta.A') <- Delta.A
  attr(out, 'Delta.T') <- Delta.T
  class(out) <- c('bbssr_powerbssr', 'data.frame')
  out
}
