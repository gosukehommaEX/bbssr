#' Simulated Rejections of a Stratified Re-Estimation Design with Random Sizes
#'
#' Internal helper simulating \code{nsim} trials of a stratified design with blinded sample
#' size re-estimation whose stratum and group sizes are random, for one scenario, and
#' \code{nsim} trials of the fixed design with the planned total and the same allocation.
#' The strata of the patients are drawn with the probabilities \code{frac}. Within stratum
#' \code{j} the groups follow permuted blocks of size \code{block} with \code{g[j]}
#' patients of group 1 (\code{alloc = 'stratum'}), or each patient is in group 1 with the
#' probability \code{b[j]} (\code{alloc = 'both'}). The interim analysis has the first
#' \code{n} patients, drawn again until every stratum has one of them, and the final
#' analysis the first \code{N} patients, where \code{N} is the total of
#' \code{strat_alloc_reestimate} for the interim stratum sizes and pooled numbers of
#' responders. The groups of the second stage continue the permuted blocks of the interim
#' patients. A group without patients is allowed in both designs, and a stratum without
#' patients in the fixed design. The random numbers come from the generator of R, so the
#' results follow \code{set.seed()}.
#'
#' @param n Interim total sample size
#' @param N.plan Planned total sample size, that of the fixed design
#' @param al List returned by \code{check_strat_alloc}
#' @param alloc \code{'stratum'} or \code{'both'}
#' @param p1 Response probabilities of group 1 in the strata
#' @param p2 Response probabilities of group 2 in the strata
#' @param nsim Number of simulated trials of each design
#' @param alpha Level of the final analysis
#' @param a List returned by \code{check_strat_args} for the final analysis
#' @param re.args List of the arguments of \code{strat_reestimate}, with the planned
#'   fractions in \code{frac}
#' @param ss.frac \code{'observed'} or \code{'planned'}
#' @param restricted Logical. Whether the final total is at least \code{N.plan}
#' @param N.min Smallest final total, or \code{NULL}
#' @param N.max Largest final total, or \code{NULL}
#'
#' @return A list with the logical vector \code{reject} and the integer vector \code{N} of
#'   final total sample sizes of the simulated trials of the design with re-estimation,
#'   and the logical vector \code{reject.trad} of the trials of the fixed design
#'
#' @keywords internal
#' @noRd
#' @importFrom stats rbinom rhyper rmultinom
strat_alloc_simulate <- function(n, N.plan, al, alloc, p1, p2, nsim, alpha, a, re.args,
                                 ss.frac, restricted, N.min, N.max) {
  J <- length(al$frac)
  B <- al$block
  g <- if (alloc == 'stratum') matrix(al$g, nsim, J, byrow = TRUE)
  # Stratum sizes of k trials of size patients, one row per trial
  strata <- function(size, k) matrix(t(rmultinom(k, size, al$frac)), k, J)
  # Size of group 1 among the first m patients of every stratum
  first <- function(m) {
    G <- if (alloc == 'stratum') {
      g * (m %/% B) + rhyper(length(m), g, B - g, m %% B)
    } else {
      rbinom(length(m), m, rep(al$b, each = nrow(m)))
    }
    matrix(as.integer(G), nrow(m), J)
  }
  # Size of group 1 among the first M patients, given G among the first m <= M. Within the
  # block of patient m, G - g floor(m / B) of its first m mod B places are of group 1; a
  # later block starts afresh
  more <- function(G, m, M) {
    if (alloc == 'both') {
      return(matrix(G + as.integer(rbinom(length(m), M - m, rep(al$b, each = nrow(m)))),
                    nrow(m), J))
    }
    out <- G
    same <- M %/% B == m %/% B
    if (any(same)) {
      h <- G[same] - g[same] * (m[same] %/% B)
      left <- B - m[same] %% B
      out[same] <- G[same] + rhyper(sum(same), g[same] - h, left - (g[same] - h),
                                    M[same] - m[same])
    }
    if (any(!same)) {
      out[!same] <- g[!same] * (M[!same] %/% B) +
        rhyper(sum(!same), g[!same], B - g[!same], M[!same] %% B)
    }
    matrix(as.integer(out), nrow(m), J)
  }
  binom <- function(size, prob) {
    matrix(as.integer(rbinom(length(size), size, rep(prob, each = nrow(size)))),
           nrow(size), J)
  }
  # Interim analysis, drawn again until every stratum has a patient
  m <- strata(n, nsim)
  bad <- which(rowSums(m == 0) > 0)
  for (try in 1:1000) {
    if (length(bad) == 0) break
    m[bad, ] <- strata(n, length(bad))
    bad <- bad[rowSums(m[bad, , drop = FALSE] == 0) > 0]
  }
  if (length(bad) > 0) {
    stop('the interim analyses of 1000 draws had no patient in one of the strata')
  }
  G <- first(m)
  A <- binom(G, p1)
  C <- binom(m - G, p2)
  S <- A + C
  # Final totals, computed once for every interim stratum sizes and pooled responders
  N <- integer(nsim)
  for (rows in split(seq_len(nsim), apply(m, 1, paste, collapse = '_'))) {
    key <- apply(S[rows, , drop = FALSE], 1, paste, collapse = '_')
    u <- !duplicated(key)
    res <- strat_alloc_reestimate(m[rows[1], ], re.args, ss.frac, N.plan, restricted,
                                  N.min, N.max, a$alternative, a$tsmethod,
                                  s = S[rows[u], , drop = FALSE])
    N[rows] <- res$N[match(key, key[u])]
  }
  # Second stage
  M <- m
  for (rows in split(seq_len(nsim), N - n)) {
    e <- N[rows[1]] - n
    if (e > 0) M[rows, ] <- M[rows, ] + strata(e, length(rows))
  }
  G2 <- more(G, m, M)
  X1 <- A + binom(G2 - G, p1)
  X2 <- C + binom((M - G2) - (m - G), p2)
  pv <- strat_pvalue_cases(G2, M - G2, X1, X2, a$type, a$alt, a$ts)
  # Fixed design with the planned total
  Mt <- strata(N.plan, nsim)
  Gt <- first(Mt)
  pt <- strat_pvalue_cases(Gt, Mt - Gt, binom(Gt, p1), binom(Mt - Gt, p2), a$type,
                           a$alt, a$ts)
  list(reject = pv %<<% alpha, N = as.integer(N), reject.trad = pt %<<% alpha)
}
