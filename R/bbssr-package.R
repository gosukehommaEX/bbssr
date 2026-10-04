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
#' @section Reuse of rejection regions:
#' A rejection region depends only on the sample sizes, the level and the test. The
#' functions that evaluate power, sample size and re-estimation designs therefore compute
#' each rejection region once and keep it for the rest of the R session, which makes a
#' second evaluation of a related design much faster than the first. The stored regions
#' are discarded when they would occupy more than about 100 MB, and
#' \code{options(bbssr.cache = FALSE)} turns the reuse off. The results are the same
#' either way.
#'
#' @keywords internal
#' @useDynLib bbssr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom utils globalVariables
"_PACKAGE"

# Column names used inside aes() are looked up in the plotted data frame rather than in
# the enclosing environment, so they are declared here to keep the check quiet
globalVariables(c('x1', 'x2', 'Reject', 'Power', 'p1', 'N2', 'Design', 'p'))

# Session store of rejection regions used by get_rr(). The limit is a number of cells of
# the outcome grid, each stored as a 4-byte logical value
.bbssr_cache <- new.env(parent = emptyenv())
.bbssr_cache$rr <- new.env(parent = emptyenv())
.bbssr_cache$cells <- 0
.bbssr_cache$max.cells <- 2.5e7
