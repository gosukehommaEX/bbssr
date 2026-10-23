#' Validate the Stratum Fractions and Allocation Ratios
#'
#' Internal helper checking the fractions of the strata and the allocation ratios of a
#' stratified design. The fractions are rescaled to sum to 1, so that they may also be
#' given as numbers of patients, and a single allocation ratio applies to every stratum.
#'
#' @param frac Fractions of the strata, positive, one element per stratum
#' @param r Allocation ratio of group 1 to group 2, a single positive value or one per
#'   stratum
#' @param J Number of strata
#'
#' @return A list with the fractions \code{frac}, summing to 1, and the vector \code{r} of
#'   length \code{J}
#'
#' @keywords internal
#' @noRd
check_strat_frac <- function(frac, r, J) {
  if (!is.numeric(frac) || length(frac) != J || anyNA(frac) || any(!is.finite(frac)) ||
      any(frac <= 0)) {
    stop('frac must be a vector of positive values with one element per stratum')
  }
  if (!is.numeric(r) || !(length(r) %in% c(1, J)) || anyNA(r) || any(!is.finite(r)) ||
      any(r <= 0)) {
    stop('r must be a single positive value or one positive value per stratum')
  }
  list(frac = frac / sum(frac), r = rep_len(as.double(r), J))
}
