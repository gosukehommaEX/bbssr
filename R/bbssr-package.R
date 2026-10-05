#' bbssr: Blinded Sample Size Re-Estimation for Binary Endpoints
#'
#' Tools for blinded sample size re-estimation (BSSR) in two-arm clinical trials with
#' binary endpoints, together with the exact power and sample size calculations the
#' re-estimation relies on. Seven tests are supported, each available with a one-sided or a
#' two-sided alternative. Two of them also test non-inferiority with a margin, and the
#' exact unconditional tests can be combined with the Berger-Boos procedure.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{BinaryRR}}}{Rejection region of an exact test}
#'   \item{\code{\link{BinaryPower}}}{Exact power at a given sample size}
#'   \item{\code{\link{BinarySampleSize}}}{Sample size attaining a target power}
#'   \item{\code{\link{BinaryPowerBSSR}}}{Power and final sample size of a BSSR design}
#'   \item{\code{\link{BinaryTypeIErrorBSSR}}}{Type I error rate of a BSSR design over the
#'     common response probability}
#'   \item{\code{\link{BinaryAlphaAdjBSSR}}}{Adjusted significance level that controls the
#'     type I error rate of a BSSR design}
#'   \item{\code{\link{BinaryBSSR}}}{Sample size re-estimation from observed interim data}
#' }
#'
#' @section Reuse of computed p-values:
#' The p-values of a test over the outcome grid depend only on the sample sizes and the
#' test, and a rejection region at any level follows from them by comparison. The
#' functions of the package therefore compute each matrix of p-values once and keep it
#' for the rest of the R session, which makes a second evaluation of a related design, or
#' a search over the significance level, much faster than the first evaluation. The stored
#' matrices are discarded when they would occupy more than about 100 MB, and
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
globalVariables(c('x1', 'x2', 'Reject', 'Power', 'p1', 'N2', 'Design', 'p', 'theta',
                  'TIE'))

# Session store of p-value matrices used by get_pvalue(). The limit is a number of cells
# of the outcome grid, each stored as an 8-byte double
.bbssr_cache <- new.env(parent = emptyenv())
.bbssr_cache$pv <- new.env(parent = emptyenv())
.bbssr_cache$cells <- 0
.bbssr_cache$max.cells <- 1.25e7
