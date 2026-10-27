#' bbssr: Blinded Sample Size Re-Estimation for Binary Endpoints
#'
#' Tools for blinded sample size re-estimation (BSSR) in two-arm clinical trials with
#' binary endpoints, together with the exact power and sample size calculations the
#' re-estimation relies on. Seven tests are supported, each available with a one-sided or a
#' two-sided alternative. Two of them also test non-inferiority with a margin on the scale
#' of the risk difference or the risk ratio, and the exact unconditional tests can be
#' combined with the Berger-Boos procedure. The same calculations are available for
#' matched-pairs designs, such as cross-over trials, with six tests based on the
#' discordant pairs, and for designs stratified into two or more strata, with the
#' Cochran-Mantel-Haenszel test, Cochran's test and the stratified exact conditional test.
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
#'   \item{\code{\link{BinaryCondRejectBSSR}}}{Rejection probabilities of a BSSR design
#'     given the pooled numbers of responders, and its type I error rate decomposed by
#'     the interim outcome}
#'   \item{\code{\link{BinaryGridBSSR}}}{Power and final sample size of several BSSR
#'     designs over common scenarios, as one data frame}
#'   \item{\code{\link{BinaryBSSR}}}{Sample size re-estimation from observed interim data}
#' }
#'
#' @section Paired binary data:
#' \describe{
#'   \item{\code{\link{PairedRR}}}{Rejection region of a test for matched pairs}
#'   \item{\code{\link{PairedPower}}}{Exact power at a given number of pairs}
#'   \item{\code{\link{PairedSampleSize}}}{Number of pairs attaining a target power}
#'   \item{\code{\link{PairedPowerBSSR}}}{Power and final number of pairs of a
#'     re-estimation design}
#'   \item{\code{\link{PairedBSSR}}}{Re-estimation from the interim number of discordant
#'     pairs}
#'   \item{\code{\link{PairedTypeIErrorBSSR}}}{Type I error rate of a re-estimation design
#'     over the probability of a discordant pair}
#'   \item{\code{\link{PairedAlphaAdjBSSR}}}{Adjusted significance level of a
#'     re-estimation design}
#' }
#'
#' @section Stratified designs:
#' \describe{
#'   \item{\code{\link{StratRR}}}{Critical values of a stratified test}
#'   \item{\code{\link{StratPower}}}{Exact power of a stratified test}
#'   \item{\code{\link{StratSampleSize}}}{Total sample size attaining a target power}
#'   \item{\code{\link{StratPowerBSSR}}}{Power and final sample size of a re-estimation
#'     design}
#'   \item{\code{\link{StratGridBSSR}}}{Power and final sample size of several
#'     re-estimation designs over common scenarios}
#'   \item{\code{\link{StratBSSR}}}{Re-estimation from the interim pooled numbers of
#'     responders of the strata}
#'   \item{\code{\link{StratTypeIErrorBSSR}}}{Type I error rate of a re-estimation
#'     design over the common response probabilities of the strata}
#'   \item{\code{\link{StratAlphaAdjBSSR}}}{Adjusted significance level of a
#'     re-estimation design with two strata}
#'   \item{\code{\link{StratCondRejectBSSR}}}{Conditional rejection probabilities of a
#'     re-estimation design with two strata, and its type I error rate decomposed by the
#'     interim outcome}
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
                  'TIE', 's', 's2', 'value', 'TIE.s', 'grp', 'Boundary', 'a', 'b', 'psi',
                  'theta1', 'theta2'))

# Session store of p-value matrices used by get_pvalue(). The limit is a number of cells
# of the outcome grid, each stored as an 8-byte double
.bbssr_cache <- new.env(parent = emptyenv())
.bbssr_cache$pv <- new.env(parent = emptyenv())
.bbssr_cache$cells <- 0
.bbssr_cache$max.cells <- 1.25e7
