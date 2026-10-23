#' Blinded Sample Size Re-Estimation for a Stratified Trial
#'
#' Re-estimates the final group sizes of the strata of an ongoing two-arm trial with a
#' binary endpoint from the blinded data of an interim analysis: the numbers of patients
#' and the pooled numbers of responders of the strata.
#'
#' @param n1 Numbers of patients of group 1 observed in the strata at the interim analysis
#' @param n2 Numbers of patients of group 2 observed in the strata at the interim analysis
#' @param S Numbers of responders observed in the strata at the interim analysis, pooled
#'   over both groups, one element per stratum
#' @param N1 Planned sizes of group 1 in the strata, or \code{NULL} (default). Required when
#'   \code{restricted} is \code{TRUE}, when \code{frac} is \code{NULL}, and when the
#'   blinded data show no effect in the direction of the alternative
#' @param N2 Planned sizes of group 2 in the strata, or \code{NULL} (default)
#' @param frac Fractions of the strata used to divide the re-estimated total sample size,
#'   see \code{\link{StratSampleSize}}. Default is \code{NULL}, which takes the fractions
#'   of the planned design
#' @inheritParams StratPowerBSSR
#'
#' @return An object of class \code{bbssr_sbssr}, a data frame with one row per stratum
#'   containing:
#' \describe{
#'   \item{Stratum}{Number of the stratum}
#'   \item{n1, n2, S}{Interim data}
#'   \item{hat.p}{Blinded estimate of the pooled response probability}
#'   \item{hat.p1, hat.p2}{Planning probabilities of the two groups, truncated to the unit
#'     interval}
#'   \item{N1.final, N2.final}{Final group sizes}
#'   \item{n1.stage2, n2.stage2}{Numbers of patients still to be enrolled}
#' }
#' The attributes \code{n.raw} (unrounded total of the normal approximation),
#' \code{N.final} (final total sample size) and \code{Power} (exact power of the final
#' design at the planning probabilities) summarize the re-estimation.
#'
#' @details
#' The re-estimation is that of \code{\link{StratPowerBSSR}}: in each stratum the blinded
#' estimate of the pooled response probability is divided between the groups with the
#' assumed effect \code{Delta.A}, and the total sample size of \code{\link{StratSampleSize}}
#' for these probabilities is divided between the strata. The final group sizes are at
#' least the interim sizes, at least the planned sizes when \code{restricted = TRUE}, and
#' their total lies between \code{N.min} and \code{N.max}.
#'
#' @examples
#' StratBSSR(n1 = c(15, 15), n2 = c(15, 15), S = c(12, 18), Delta.A = 0.2, r = 1,
#'           alpha = 0.025, tar.power = 0.8, Test = 'CMH', N1 = c(30, 30),
#'           N2 = c(30, 30), ss.method = 'standard')
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @seealso \code{\link{StratPowerBSSR}}
#' @export
StratBSSR <- function(n1, n2, S, Delta.A, r, alpha, tar.power, Test, restricted = FALSE,
                      N1 = NULL, N2 = NULL, frac = NULL,
                      alternative = c('greater', 'less', 'two.sided'),
                      tsmethod = c('minlike', 'central', 'blaker'),
                      effect = c('RD', 'RR', 'OR'),
                      ss.method = c('exact', 'standard'),
                      ss.Test = Test, ss.alpha = alpha,
                      search = c('crossing', 'smallest', 'stable'),
                      search.limit = c(2, 50), N.min = NULL, N.max = NULL) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  effect <- match.arg(effect)
  ss.method <- match.arg(ss.method)
  search <- match.arg(search)
  a <- check_strat_args(n1, n2, alpha, Test, alternative, tsmethod)
  J <- length(a$N1)
  if (!is.numeric(S) || length(S) != J || anyNA(S) || any(S != round(S)) || any(S < 0) ||
      any(S > a$N1 + a$N2)) {
    stop('S must have one integer per stratum between 0 and n1 + n2')
  }
  if (!is.null(N1) || !is.null(N2)) {
    plan <- check_strat_args(N1, N2, alpha, Test, alternative, tsmethod)
    if (length(plan$N1) != J) stop('N1 and N2 must have one element per stratum')
    N1 <- plan$N1
    N2 <- plan$N2
  }
  if (is.null(frac)) {
    if (is.null(N1)) stop('supply frac, or the planned sizes N1 and N2')
    frac <- N1 + N2
  }
  fr <- check_strat_frac(frac, r, J)
  re <- strat_reestimate(matrix(S, nrow = 1), a$N1, a$N2, Delta.A, effect, fr$frac, fr$r,
                         ss.method, ss.Test, ss.alpha, tar.power, a$alternative,
                         a$tsmethod, N1, N2, restricted, N.min, N.max, search,
                         search.limit)
  N1.final <- as.integer(re$N1[1, ])
  N2.final <- as.integer(re$N2[1, ])
  Power <- strat_fixed_power(N1.final, N2.final, re$hat.p1, re$hat.p2, alpha, a$Test,
                             a$alternative, a$tsmethod)
  out <- data.frame(Stratum = seq_len(J), n1 = a$N1, n2 = a$N2, S = as.integer(S),
                    hat.p = re$hat.p[1, ], hat.p1 = re$hat.p1[1, ],
                    hat.p2 = re$hat.p2[1, ], N1.final = N1.final, N2.final = N2.final,
                    n1.stage2 = N1.final - a$N1, n2.stage2 = N2.final - a$N2)
  attr(out, 'n.raw') <- re$n.raw
  attr(out, 'N.final') <- as.integer(sum(N1.final) + sum(N2.final))
  attr(out, 'Power') <- Power
  attr(out, 'Test') <- a$Test
  attr(out, 'alternative') <- a$alternative
  attr(out, 'tsmethod') <- a$tsmethod
  attr(out, 'alpha') <- alpha
  attr(out, 'tar.power') <- tar.power
  attr(out, 'restricted') <- restricted
  attr(out, 'Delta.A') <- Delta.A
  attr(out, 'effect') <- effect
  attr(out, 'frac') <- fr$frac
  attr(out, 'r') <- fr$r
  attr(out, 'ss.method') <- ss.method
  attr(out, 'ss.Test') <- ss.Test
  attr(out, 'ss.alpha') <- ss.alpha
  attr(out, 'search') <- search
  attr(out, 'search.limit') <- re$limit
  attr(out, 'N.min') <- N.min
  attr(out, 'N.max') <- N.max
  class(out) <- c('bbssr_sbssr', 'data.frame')
  out
}
