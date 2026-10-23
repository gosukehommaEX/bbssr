# Reference computations for the paired design, by direct enumeration of the multinomial
# outcomes. They are slow and are used with small numbers of pairs only.

# Probabilities of the outcomes (a, b) of N pairs, with a pairs in which only treatment 1
# responds and b pairs in which only treatment 2 responds, as an (N + 1) by (N + 1) matrix
trinomial_ref <- function(N, p10, p01) {
  P <- matrix(0, N + 1, N + 1)
  for (a in 0:N) {
    for (b in 0:(N - a)) {
      P[a + 1, b + 1] <- stats::dmultinom(c(a, b, N - a - b),
                                          prob = c(p10, p01, max(0, 1 - p10 - p01)))
    }
  }
  P
}

# Rejection probability of a re-estimation design by enumeration of the interim and
# second-stage outcomes. sizes[D1 + 1] is the final number of pairs after D1 interim
# discordant pairs and regions a list of regions named by the number of pairs
paired_bssr_ref <- function(sizes, regions, n1, p10, p01) {
  T1 <- trinomial_ref(n1, p10, p01)
  tot <- 0
  for (a1 in 0:n1) {
    for (b1 in 0:(n1 - a1)) {
      N <- sizes[a1 + b1 + 1]
      R <- regions[[as.character(N)]]
      T2 <- trinomial_ref(N - n1, p10, p01)
      sub <- R[a1 + 1:(N - n1 + 1), b1 + 1:(N - n1 + 1)]
      tot <- tot + T1[a1 + 1, b1 + 1] * sum(T2 * sub)
    }
  }
  tot
}

# A random rejection region of N pairs, FALSE at the impossible outcomes
random_region <- function(N, prob = 0.5) {
  R <- matrix(stats::runif((N + 1)^2) < prob, N + 1, N + 1)
  R[outer(0:N, 0:N, '+') > N] <- FALSE
  R
}
