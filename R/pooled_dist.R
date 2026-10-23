#' Distribution of the Pooled Number of Responders of Two Groups
#'
#' Internal helper returning the distribution of the sum of independent binomial numbers
#' of responders with \code{n1} trials and probability \code{p1}, and \code{n2} trials and
#' probability \code{p2}.
#'
#' @param n1 Size of group 1
#' @param n2 Size of group 2
#' @param p1 Response probability of group 1
#' @param p2 Response probability of group 2
#'
#' @return A numeric vector of the probabilities of \code{0, ..., n1 + n2}
#'
#' @keywords internal
#' @noRd
#' @importFrom stats dbinom
pooled_dist <- function(n1, n2, p1, p2) {
  w <- outer(dbinom(0:n1, n1, p1), dbinom(0:n2, n2, p2))
  as.vector(rowsum(as.vector(w), as.vector(outer(0:n1, 0:n2, '+')), reorder = TRUE))
}
