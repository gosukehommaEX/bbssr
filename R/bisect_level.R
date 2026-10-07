#' Bisection over the Nominal Level with Certification of the Result
#'
#' Internal helper of \code{BinaryAlphaAdjBSSR} finding the largest nominal level at which
#' the largest type I error rate of a design does not exceed \code{alpha}, for designs
#' whose rejection regions shrink as the level decreases, so that the largest rate is a
#' non-decreasing function of the level.
#'
#' A level passes when the largest rate \code{y} is at most \code{alpha} and, when an
#' upper bound \code{bound} of the largest rate is available, so is the bound. The nominal
#' level \code{alpha} is assessed and, with \code{certify}, certified first. If it fails,
#' the interval from 0 to \code{alpha} is bisected until its width is at most
#' \code{tol * alpha}, with each level assessed by \code{assess}. The value of the
#' response probability at which the last failing level exceeded \code{alpha} is passed
#' to \code{assess} as a probe. With \code{certify}, the level found is then certified.
#' If the certified maximum or its upper bound exceeds \code{alpha}, the assessment has
#' missed a value of the response probability at which the level fails, or the maximum
#' lies within the tolerance of the bound below \code{alpha}, and the bisection is
#' repeated between 0 and that level with every level that passes the assessment also
#' certified, so the result is certified in either case.
#'
#' @param make Function of the level returning a design, which is passed to
#'   \code{assess} and \code{certify}
#' @param assess Function of a design and a probe (a value of the response probability,
#'   or \code{NULL}) returning a list with the location \code{x}, the largest rate
#'   \code{y} and \code{bound = NA}
#' @param certify \code{NULL}, or a function of a design and the result of \code{assess}
#'   returning that result with the certified maximum and its upper bound \code{bound}
#' @param alpha Target level
#' @param tol Tolerance of the bisection, relative to \code{alpha}
#'
#' @return A list with the level \code{level} found, the result \code{m0} at the nominal
#'   level and the result \code{m} at the level found
#'
#' @keywords internal
#' @noRd
bisect_level <- function(make, assess, certify, alpha, tol) {
  passes <- function(m) m$y <= alpha && (is.na(m$bound) || m$bound <= alpha)
  des <- make(alpha)
  m0 <- assess(des, NULL)
  if (!is.null(certify)) m0 <- certify(des, m0)
  if (passes(m0)) return(list(level = alpha, m0 = m0, m = m0))
  # Result for the level 0, at which nothing is rejected
  none <- list(x = NA_real_, y = 0, bound = if (is.null(certify)) NA_real_ else 0)
  probe <- m0$x
  lo <- 0
  hi <- alpha
  m.lo <- none
  strict <- FALSE
  repeat {
    while (hi - lo > tol * alpha) {
      mid <- (lo + hi) / 2
      des <- make(mid)
      m <- assess(des, probe)
      if (strict && m$y <= alpha) m <- certify(des, m)
      if (passes(m)) {
        lo <- mid
        m.lo <- m
      } else {
        hi <- mid
        probe <- m$x
      }
    }
    if (is.null(certify) || strict || lo == 0) break
    m.lo <- certify(make(lo), m.lo)
    if (passes(m.lo)) break
    # The level fails after all. Bisect again below it, certifying every level that passes
    hi <- lo
    lo <- 0
    probe <- m.lo$x
    m.lo <- none
    strict <- TRUE
  }
  list(level = lo, m0 = m0, m = m.lo)
}
