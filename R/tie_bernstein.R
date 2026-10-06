#' Bernstein Coefficients of the Type I Error Rate over an Interval
#'
#' Internal helper returning the coefficients of the rejection probability of a design in
#' the Bernstein basis \code{B_{k,n}(t) = choose(n, k) t^k (1 - t)^(n - k)} of the
#' variable \code{t} in \code{[0, 1]}, when the response probabilities run linearly from
#' \code{p1[1]} and \code{p2[1]} at \code{t = 0} to \code{p1[2]} and \code{p2[2]} at
#' \code{t = 1}.
#'
#' A binomial probability in a response probability \code{a (1 - t) + c t} is
#' \code{b(X; N, a (1 - t) + c t) = sum_j M[X, j] B_{j,N}(t)}, where column \code{j} of
#' \code{M}, computed by \code{binom_conv_matrix}, is the distribution of the sum of
#' independent \code{Bin(N - j, a)} and \code{Bin(j, c)} variables. On \code{[0, 1]}, with
#' \code{a = 0} and \code{c = 1}, \code{M} is the identity. The product
#' \code{B_{j,N1}(t) B_{l,N2}(t)} equals \code{dhyper(j, N1, N2, j + l) B_{j+l,N1+N2}(t)},
#' so each final sample size contributes a polynomial of degree \code{N1 + N2}. The
#' contributions are added in increasing order of the degree, which is raised one step at a
#' time by \code{c'_i = (i / (n + 1)) c_{i-1} + (1 - i / (n + 1)) c_i}. Every quantity in
#' these steps is non-negative, so the coefficients are free of cancellation.
#'
#' @param weights List returned by \code{tie_weights}
#' @param p1 Response probabilities of group 1 at the two ends of the interval
#' @param p2 Response probabilities of group 2 at the two ends of the interval
#'
#' @return The coefficients, a numeric vector of length one more than the largest
#'   \code{N1 + N2}
#'
#' @keywords internal
#' @noRd
#' @importFrom stats dhyper
tie_bernstein <- function(weights, p1, p2) {
  deg <- vapply(weights, function(w) as.numeric(w$N1 + w$N2), numeric(1))
  coef <- 0
  for (k in order(deg)) {
    w <- weights[[k]]
    D <- w$C
    if (p1[1] != 0 || p1[2] != 1) D <- crossprod(binom_conv_matrix(w$N1, p1[1], p1[2]), D)
    if (p2[1] != 0 || p2[2] != 1) D <- D %*% binom_conv_matrix(w$N2, p2[1], p2[2])
    S <- outer(0:w$N1, 0:w$N2, '+')
    h <- rowsum(as.vector(D) * dhyper(as.vector(row(D)) - 1, w$N1, w$N2, as.vector(S)),
                as.vector(S))
    n <- length(coef) - 1
    while (n < deg[k]) {
      i <- 0:(n + 1)
      coef <- c(coef, 0) * (1 - i / (n + 1)) + c(0, coef) * (i / (n + 1))
      n <- n + 1
    }
    coef <- coef + as.vector(h)
  }
  coef
}
