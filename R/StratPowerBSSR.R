#' Power of Blinded Sample Size Re-Estimation for a Stratified Design
#'
#' Calculates the power and the distribution of the final total sample size of a two-arm
#' trial with a binary endpoint, stratified into two or more strata, with blinded sample
#' size re-estimation, together with the power of the fixed design with the planned group
#' sizes. The sample size is re-estimated at an interim analysis from the pooled numbers
#' of responders of the strata, without the treatment groups, and the final analysis is a
#' test of \code{\link{StratRR}} with all patients. With two strata the power is computed
#' exactly, and with three or more strata by simulation.
#'
#' @param p True pooled response probabilities of the strata: a vector with one element
#'   per stratum, or a matrix with one row per scenario and one column per stratum
#' @param Delta.A Assumed treatment effect, on the scale given by \code{effect}, common to
#'   the strata and used to divide the blinded estimate of the pooled response probability
#'   of each stratum between the groups at the re-estimation
#' @param Delta.T True treatment effect, on the scale given by \code{effect}: a single value
#'   common to the strata, or one value per stratum
#' @param N1 Planned sizes of group 1 in the strata, which are those of the fixed design
#' @param N2 Planned sizes of group 2 in the strata
#' @param omega Fraction of the planned sizes at the interim analysis: in stratum \code{j}
#'   the interim size of group 2 is \code{ceiling(omega N2[j])} and that of group 1 is
#'   \code{ceiling(r[j] n2[j])}. Supply either \code{omega} or \code{n.interim}
#' @param r Allocation ratio of group 1 to group 2 within the strata, a single value or one
#'   value per stratum. The blinded estimate of the pooled response probability of a
#'   stratum is divided between the groups with this ratio, also when the interim sizes
#'   are not in this ratio
#' @param alpha Level of significance of the final analysis
#' @param tar.power Target power of the re-estimation
#' @param restricted Logical. If \code{TRUE}, the final group sizes are at least the planned
#'   sizes. Default is \code{FALSE}
#' @param effect Scale of \code{Delta.A} and \code{Delta.T}. Options: \code{'RD'} (default)
#'   for the risk difference, \code{'RR'} for the risk ratio or \code{'OR'} for the odds
#'   ratio
#' @param ss.method Re-estimation rule. \code{'exact'} (default) re-estimates the sample
#'   size of \code{\link{StratSampleSize}} with \code{method = 'exact'}, and
#'   \code{'standard'} that of the normal approximation. With three or more strata only
#'   \code{'standard'} is available
#' @param ss.Test Test used in the re-estimation. Default is \code{Test}
#' @param ss.alpha Level used in the re-estimation. Default is \code{alpha}
#' @param ss.strata Blinded estimate of the re-estimation. \code{'stratified'} (default)
#'   estimates the pooled response probability of each stratum from its own interim
#'   patients, and \code{'pooled'} uses the proportion of responders among all the interim
#'   patients for every stratum, as a re-estimation that ignores the strata would. With
#'   \code{'pooled'} the final group sizes depend only on the total number of interim
#'   responders
#' @param search Search of the exact re-estimation, see \code{\link{StratSampleSize}}
#' @param search.limit Limit of the search \code{'stable'}, see
#'   \code{\link{StratSampleSize}}
#' @param N.min Smallest final total sample size, or \code{NULL} (default)
#' @param N.max Largest final total sample size, or \code{NULL} (default)
#' @param n.interim Interim group sizes, as a matrix with two rows (group 1 and group 2)
#'   and one column per stratum. Supply either \code{omega} or \code{n.interim}
#' @param frac Fractions of the strata used to divide the re-estimated total sample size
#'   between the strata, see \code{\link{StratSampleSize}}. Default is \code{NULL}, which
#'   takes the fractions of the planned design, \code{(N1 + N2) / sum(N1 + N2)}
#' @param nsim Number of simulated trials for each scenario with three or more strata.
#'   Default is \code{1e5}. Ignored with two strata
#' @param ss.frac Fractions of the strata used by the re-estimation with a random
#'   allocation: \code{'observed'} (default), those of the interim patients, or
#'   \code{'planned'}, those given by \code{frac}
#' @param nseq Number of allocation sequences drawn with a random allocation and two
#'   strata. Default is 100
#' @param alloc Allocation of the patients to the strata and the groups. \code{'fixed'}
#'   (default) uses the planned and interim group sizes. The others make the sizes random,
#'   see Details: \code{'stratum'} (stratum sizes random, permuted blocks within the
#'   strata) and \code{'both'} (stratum and group sizes random), with any number of
#'   strata, and \code{'group'} (stratum sizes in the fixed fractions, group sizes random),
#'   with two strata
#' @inheritParams StratPower
#' @inheritParams StratRR
#'
#' @return An object of class \code{bbssr_spowerbssr}, a data frame with one row for each
#'   scenario in which both response probabilities of every stratum lie in the unit
#'   interval, containing:
#' \describe{
#'   \item{p.S1, p.S2, ...}{True pooled response probabilities of the strata}
#'   \item{p1.S1, ..., p2.S1, ...}{True response probabilities of the two groups}
#'   \item{power.BSSR}{Power of the design with re-estimation}
#'   \item{power.TRAD}{Power of the fixed design with the planned group sizes}
#'   \item{E.N}{Expected final total sample size}
#'   \item{se.BSSR}{Monte Carlo standard error of \code{power.BSSR}, with three or more
#'     strata or with a random allocation}
#'   \item{se.TRAD}{Monte Carlo standard error of \code{power.TRAD}, with three or more
#'     strata and a random allocation}
#' }
#' The attribute \code{reestimation} holds, with two strata, the final group sizes for every
#' pair of interim pooled numbers of responders, \code{N.dist} the distribution of the
#' final total sample size in each scenario, and the other attributes the design settings.
#'
#' @details
#' In stratum \code{j} the blinded estimate of the pooled response probability is the
#' number \code{s_j} of interim responders divided by the interim size of the stratum. It
#' is divided between the groups with the assumed effect \code{Delta.A} and the allocation
#' ratio of the stratum, as in \code{\link{BinaryPowerBSSR}}, and the total sample size is
#' re-estimated from these probabilities by the rule \code{ss.method}. The total is divided
#' between the strata and the groups as in \code{\link{StratSampleSize}}, and the final
#' group sizes are at least the interim sizes, at least the planned sizes when
#' \code{restricted = TRUE}, with a total of at least \code{N.min} and at most
#' \code{N.max}. When the probabilities of every stratum show no effect in the direction of
#' the alternative, which happens with \code{effect = 'RR'} when no interim patient
#' responds and with \code{effect = 'OR'} when in each stratum no interim patient or every
#' interim patient responds, the planned group sizes are kept. A stratum whose pooled
#' estimate is close to 0 with \code{effect = 'RR'}, or close to 0 or 1 with
#' \code{effect = 'OR'}, has planning probabilities close to each other, so the
#' re-estimated size can be very large, and the time of the exact computation grows with
#' it; an upper bound \code{N.max} keeps both in check.
#'
#' With two strata the final group sizes depend on the interim data only through
#' \code{(s_1, s_2)}. Given \code{s_j}, the number of interim responders of group 1 in
#' stratum \code{j} follows from the two binomial distributions of the stratum, and the
#' second-stage numbers of responders are binomial given the final sizes. The power is the
#' sum over \code{(s_1, s_2)} of the probability of the interim outcome times the
#' conditional probability of rejection, computed exactly by enumeration in compiled code.
#'
#' With three or more strata the number of interim outcomes grows as a product over the
#' strata, and the power is estimated from \code{nsim} simulated trials, with the random
#' numbers of R, so that \code{set.seed()} makes the results reproducible. The power of the
#' fixed design with the fixed allocation is exact for any number of strata.
#'
#' With a random allocation (\code{alloc} other than \code{'fixed'}) and two strata the
#' patients enter the strata and the groups as in \code{\link{StratPower}}, the interim
#' analysis has the first \code{sum(n1 + n2)} patients of the interim sizes above, and the
#' final analysis the first \code{N} patients, where \code{N} is the re-estimated total.
#' The re-estimation computes the total as above, for the fixed allocation, with the
#' fractions of the strata among the interim patients (\code{ss.frac = 'observed'}) or
#' those of \code{frac}, and takes the total of the group sizes it gives, within the same
#' bounds; the allocation of the remaining patients is random. Given an allocation
#' sequence, that is the strata and the groups of the patients in the order of enrolment,
#' the design has fixed group sizes and its power is computed exactly as above. The power
#' is the average over \code{nseq} sequences drawn with the random numbers of R, so that
#' \code{set.seed()} makes the results reproducible, and \code{se.BSSR} is the standard
#' deviation of the powers of the sequences divided by \code{sqrt(nseq)}. The power of the
#' fixed design with the planned total and the same allocation is exact. A sequence
#' without interim patients in a stratum is drawn again.
#'
#' With a random allocation and three or more strata (\code{alloc = 'stratum'} or
#' \code{'both'}) the allocation is simulated with the responses, in \code{nsim} trials.
#' The strata of the first \code{sum(n1 + n2)} patients are multinomial with the
#' probabilities \code{frac}, drawn again until every stratum has one of them, and their
#' groups follow the allocation within the strata. The total is re-estimated as with two
#' strata, and the strata and the groups of the remaining patients continue the
#' allocation, the permuted blocks of \code{alloc = 'stratum'} included. The power of the
#' fixed design with the planned total and the same allocation is estimated from
#' \code{nsim} other trials, with the standard error \code{se.TRAD}. A group without
#' patients is allowed in both designs, and the fixed design is analysed without a
#' stratum that has no patients.
#'
#' @references
#' Woolson RF, Bean JA, Rojas PB (1986). Sample size for case-control studies using
#' Cochran's statistic. \emph{Biometrics}, 42(4), 927-932.
#'
#' Nam J (1992). Sample size determination for case-control studies and the comparison of
#' stratified and unstratified analyses. \emph{Biometrics}, 48(2), 389-395.
#'
#' Jung SH (2014). Stratified Fisher's exact test and its sample size calculation.
#' \emph{Biometrical Journal}, 56(1), 129-140.
#'
#' @examples
#' res <- StratPowerBSSR(p = rbind(c(0.3, 0.5), c(0.4, 0.6)), Delta.A = 0.2,
#'                       Delta.T = 0.2, N1 = c(30, 30), N2 = c(30, 30), omega = 0.5,
#'                       r = 1, alpha = 0.025, tar.power = 0.8, Test = 'CMH',
#'                       ss.method = 'standard')
#' print(res)
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @seealso \code{\link{StratTypeIErrorBSSR}}, \code{\link{StratBSSR}}
#' @export
StratPowerBSSR <- function(p, Delta.A, Delta.T, N1, N2, omega = NULL, r, alpha, tar.power,
                           Test, restricted = FALSE,
                           alternative = c('greater', 'less', 'two.sided'),
                           tsmethod = c('minlike', 'central', 'blaker'),
                           effect = c('RD', 'RR', 'OR'),
                           ss.method = c('exact', 'standard'),
                           ss.Test = Test, ss.alpha = alpha,
                           ss.strata = c('stratified', 'pooled'),
                           search = c('crossing', 'smallest', 'stable'),
                           search.limit = c(2, 50), N.min = NULL, N.max = NULL,
                           n.interim = NULL, frac = NULL, nsim = 1e5,
                           alloc = c('fixed', 'stratum', 'group', 'both'), block = 4,
                           ss.frac = c('observed', 'planned'), nseq = 100) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  effect <- match.arg(effect)
  ss.method <- match.arg(ss.method)
  search <- match.arg(search)
  ss.strata <- match.arg(ss.strata)
  alloc <- match.arg(alloc)
  ss.frac <- match.arg(ss.frac)
  d <- strat_design(N1, N2, r, frac, omega, n.interim, alpha, Test, alternative, tsmethod,
                    ss.method, nsim)
  a <- d$a
  J <- d$J
  re.args <- list(Delta.A = Delta.A, effect = effect, frac = d$frac, r = d$r,
                  ss.method = ss.method, ss.Test = ss.Test, ss.alpha = ss.alpha,
                  tar.power = tar.power, N1.plan = a$N1, N2.plan = a$N2,
                  restricted = restricted, N.min = N.min, N.max = N.max, search = search,
                  search.limit = search.limit, ss.strata = ss.strata)
  # True probabilities of the groups. Scenarios with a probability outside the unit
  # interval are omitted, and probabilities outside it only by rounding error are moved
  # onto it
  P <- check_strat_probs(p, J, 'p')
  if (!is.numeric(Delta.T) || !(length(Delta.T) %in% c(1, J)) || anyNA(Delta.T)) {
    stop('Delta.T must be a single value or one value per stratum')
  }
  dT <- rep_len(Delta.T, J)
  p1 <- p2 <- matrix(NA_real_, nrow(P), J)
  for (j in seq_len(J)) {
    sp <- split_pooled(P[, j], dT[j], d$r[j], effect)
    p1[, j] <- sp$p1
    p2[, j] <- sp$p2
  }
  out.of.range <- (p1 %<<% 0) | (p1 %>>% 1) | (p2 %<<% 0) | (p2 %>>% 1)
  keep <- rowSums(matrix(out.of.range, ncol = J)) == 0
  if (!any(keep)) {
    stop('no scenario has the response probabilities of both groups inside the unit ',
         'interval in every stratum')
  }
  P <- P[keep, , drop = FALSE]
  p1 <- pmin(pmax(p1[keep, , drop = FALSE], 0), 1)
  p2 <- pmin(pmax(p2[keep, , drop = FALSE], 0), 1)
  K <- nrow(P)
  se.BSSR <- NULL
  se.TRAD <- NULL
  power.TRAD <- NULL
  map <- NULL
  if (alloc != 'fixed' && J > 2) {
    # The allocation is simulated with the responses
    st <- list(al = check_strat_alloc(alloc, NULL, d$frac, d$r, block, J),
               N.plan = sum(a$N1 + a$N2), n = sum(d$n1 + d$n2))
    if (!is.null(N.max) && N.max < st$n) {
      stop('N.max must be at least the interim total sample size')
    }
  }
  if (alloc != 'fixed' && J == 2) {
    st <- strat_alloc_setup(d, alloc, block, nseq, function(level) re.args, ss.frac,
                            restricted, N.min, N.max)
    lay <- st$layout(ss.alpha)
    pw <- strat_seq_power(lay$inputs, p1, p2, a$type, a$alt, a$ts, alpha)
    power.BSSR <- colMeans(pw)
    se.BSSR <- apply(pw, 2, stats::sd) / sqrt(nseq)
    # Probabilities of the interim cells of each sequence, s1 varying fastest
    prob <- function(s, k) {
      as.vector(outer(pooled_dist(s$n1[1], s$n2[1], p1[k, 1], p2[k, 1]),
                      pooled_dist(s$n1[2], s$n2[2], p1[k, 2], p2[k, 2])))
    }
    E.N <- vapply(seq_len(K), function(k) {
      mean(vapply(lay$inputs, function(s) sum(prob(s, k) * s$Nfin), numeric(1)))
    }, numeric(1))
    N.dist <- do.call(rbind, lapply(seq_len(K), function(k) {
      w <- unlist(lapply(lay$inputs, function(s) prob(s, k))) / nseq
      N <- unlist(lapply(lay$inputs, function(s) s$Nfin))
      w <- rowsum(w, N, reorder = TRUE)[, 1]
      data.frame(scenario = k, N = as.integer(names(w)), prob = unname(w))
    }))
  } else if (J == 2) {
    sm <- strat_map(d$n1, d$n2, alpha, a$Test, a$alternative, a$tsmethod, re.args)
    map <- sm$map
    power.BSSR <- strat_bssr_power(d$n1, d$n2, sm$designs, sm$id, p1, p2)
    prob <- lapply(seq_len(K), function(k) {
      as.vector(outer(pooled_dist(d$n1[1], d$n2[1], p1[k, 1], p2[k, 1]),
                      pooled_dist(d$n1[2], d$n2[2], p1[k, 2], p2[k, 2])))
    })
    E.N <- vapply(prob, function(w) sum(w * map$N), numeric(1))
    N.dist <- do.call(rbind, lapply(seq_len(K), function(k) {
      w <- rowsum(prob[[k]], map$N, reorder = TRUE)[, 1]
      data.frame(scenario = k, N = as.integer(names(w)), prob = unname(w))
    }))
  } else {
    sims <- lapply(seq_len(K), function(k) {
      if (alloc == 'fixed') {
        strat_simulate(d$n1, d$n2, p1[k, ], p2[k, ], d$nsim, alpha, a, re.args)
      } else {
        strat_alloc_simulate(st$n, st$N.plan, st$al, alloc, p1[k, ], p2[k, ], d$nsim,
                             alpha, a, re.args, ss.frac, restricted, N.min, N.max)
      }
    })
    power.BSSR <- vapply(sims, function(s) mean(s$reject), numeric(1))
    se.BSSR <- sqrt(power.BSSR * (1 - power.BSSR) / d$nsim)
    if (alloc != 'fixed') {
      power.TRAD <- vapply(sims, function(s) mean(s$reject.trad), numeric(1))
      se.TRAD <- sqrt(power.TRAD * (1 - power.TRAD) / d$nsim)
    }
    E.N <- vapply(sims, function(s) mean(s$N), numeric(1))
    N.dist <- do.call(rbind, lapply(seq_len(K), function(k) {
      w <- table(sims[[k]]$N) / d$nsim
      data.frame(scenario = k, N = as.integer(names(w)), prob = as.vector(w))
    }))
  }
  rownames(N.dist) <- NULL
  if (is.null(power.TRAD)) {
    power.TRAD <- if (alloc == 'fixed') {
      strat_fixed_power(a$N1, a$N2, p1, p2, alpha, a$Test, a$alternative, a$tsmethod)
    } else {
      strat_alloc_power(st$N.plan, st$al, alloc, p1, p2, alpha, a$Test, a$alternative,
                        a$tsmethod)$power
    }
  }
  out <- data.frame(P, p1, p2, power.BSSR = power.BSSR, power.TRAD = power.TRAD,
                    E.N = E.N)
  names(out)[seq_len(3 * J)] <- c(paste0('p.S', seq_len(J)), paste0('p1.S', seq_len(J)),
                                  paste0('p2.S', seq_len(J)))
  if (!is.null(se.BSSR)) out$se.BSSR <- se.BSSR
  if (!is.null(se.TRAD)) out$se.TRAD <- se.TRAD
  attr(out, 'Test') <- a$Test
  attr(out, 'alternative') <- a$alternative
  attr(out, 'tsmethod') <- a$tsmethod
  attr(out, 'alpha') <- alpha
  attr(out, 'tar.power') <- tar.power
  attr(out, 'restricted') <- restricted
  attr(out, 'N1') <- a$N1
  attr(out, 'N2') <- a$N2
  attr(out, 'n1.interim') <- d$n1
  attr(out, 'n2.interim') <- d$n2
  attr(out, 'frac') <- d$frac
  attr(out, 'r') <- d$r
  attr(out, 'Delta.A') <- Delta.A
  attr(out, 'Delta.T') <- Delta.T
  attr(out, 'effect') <- effect
  attr(out, 'ss.method') <- ss.method
  attr(out, 'ss.Test') <- ss.Test
  attr(out, 'ss.alpha') <- ss.alpha
  attr(out, 'ss.strata') <- ss.strata
  attr(out, 'search') <- search
  attr(out, 'N.min') <- N.min
  attr(out, 'N.max') <- N.max
  attr(out, 'nsim') <- if (J > 2) d$nsim else NULL
  attr(out, 'reestimation') <- map
  attr(out, 'N.dist') <- N.dist
  attr(out, 'alloc') <- alloc
  if (alloc != 'fixed') {
    attr(out, 'block') <- st$al$block
    attr(out, 'ss.frac') <- ss.frac
    if (J == 2) attr(out, 'nseq') <- as.integer(nseq)
    attr(out, 'n.interim') <- as.integer(st$n)
  }
  class(out) <- c('bbssr_spowerbssr', 'data.frame')
  out
}
