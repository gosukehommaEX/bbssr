#' bbssr: Blinded Sample Size Re-Estimation for Binary Endpoints
#'
#' Tools for blinded sample size re-estimation (BSSR) in two-arm clinical trials with
#' binary endpoints, together with the exact power and sample size calculations the
#' re-estimation relies on. Five exact tests are supported, each available with a one-sided
#' or a two-sided alternative, and the exact unconditional tests can be combined with the
#' Berger-Boos procedure.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{BinaryRR}}}{Rejection region of an exact test}
#'   \item{\code{\link{BinaryPower}}}{Exact power at a given sample size}
#'   \item{\code{\link{BinarySampleSize}}}{Sample size attaining a target power}
#'   \item{\code{\link{BinaryPowerBSSR}}}{Operating characteristics of a BSSR design}
#'   \item{\code{\link{BinaryBSSR}}}{Sample size re-estimation from observed interim data}
#' }
#'
#' @keywords internal
#' @useDynLib bbssr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom utils globalVariables
"_PACKAGE"

# Column names used inside aes() are looked up in the plotted data frame rather than in
# the enclosing environment, so they are declared here to keep the check quiet
globalVariables(c('x1', 'x2', 'Reject', 'Power', 'p1', 'N2', 'Design', 'p'))
