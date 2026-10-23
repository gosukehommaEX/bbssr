#' Weights of the Discordant Pairs for the Certified Type I Error Rate
#'
#' Internal helper returning, for every distinct final number of pairs \code{N_k} of a
#' design with blinded sample size re-estimation for paired binary data, the weight
#' \code{W_k(D)} of \code{D} discordant pairs among the \code{N_k} pairs. The interim
#' analysis observes \code{D1} discordant pairs among its \code{n1} pairs, and the final
#' number of pairs is \code{N_k} for the values of \code{D1} in a set \code{S_k}. Since
#' \code{b(D1; n1, psi) b(D2; N_k - n1, psi) = b(D; N_k, psi) h(D1; D)} with
#' \code{D = D1 + D2} and \code{h} the hypergeometric probability of \code{D1} discordant
#' pairs among \code{n1} pairs drawn from \code{N_k} pairs of which \code{D} are
#' discordant, the weight is \code{W_k(D) = sum over D1 in S_k of h(D1; D)}. The rejection
#' probability of the design is then the sum over \code{k} and over the outcomes
#' \code{(a, b)} of \code{R_k(a, b) W_k(a + b)} times the multinomial probability of
#' \code{(a, b)} among \code{N_k} pairs.
#'
#' @param N.D1 Final number of pairs for \code{D1 = 0, ..., n1}
#' @param n1 Number of pairs at the interim analysis
#' @param Ns Distinct final numbers of pairs
#'
#' @return A list with one numeric vector of length \code{N_k + 1} for each element of
#'   \code{Ns}
#'
#' @keywords internal
#' @noRd
#' @importFrom stats dhyper
paired_tie_weights <- function(N.D1, n1, Ns) {
  lapply(Ns, function(Nk) {
    D1 <- which(N.D1 == Nk) - 1L
    h <- outer(0:Nk, D1, function(D, d) dhyper(d, D, Nk - D, n1))
    as.vector(rowSums(h))
  })
}
