# Design with two final sample sizes and interim cells assigned to them alternately
make_design <- function(prob, seed) {
  set.seed(seed)
  n11 <- 4L
  n12 <- 3L
  n21 <- c(3L, 6L)
  n22 <- c(2L, 5L)
  rr.list <- lapply(1:2, function(k) {
    nr <- n11 + n21[k] + 1L
    nc <- n12 + n22[k] + 1L
    matrix(stats::runif(nr * nc) < prob, nrow = nr, ncol = nc)
  })
  x11 <- rep(0:n11, times = n12 + 1L)
  x12 <- rep(0:n12, each = n11 + 1L)
  list(rr.list = rr.list, rr.id = as.integer((x11 + x12) %% 2L), x11 = x11, x12 = x12,
       n21 = n21, n22 = n22, n11 = n11, n12 = n12)
}

test_that("bssr_power reproduces a direct summation for arbitrary rejection regions", {
  d <- make_design(0.4, 1)
  # The regions must contain columns with several runs, or the test would only cover
  # the one-sided case
  expect_true(any(unlist(lapply(d$rr.list, column_run_count)) >= 2))
  p1 <- c(0.2, 0.55, 0.9)
  p2 <- c(0.1, 0.5, 0.3)
  got <- bssr_power(d$rr.list, d$rr.id, d$x11, d$x12, d$n21, d$n22, p1, p2,
                    d$n11, d$n12)
  want <- bssr_power_ref(d$rr.list, d$rr.id, d$x11, d$x12, d$n21, d$n22, p1, p2,
                         d$n11, d$n12)
  expect_equal(got, want, tolerance = 1e-12)
})

test_that("bssr_power gives one for a full region and zero for an empty one", {
  d <- make_design(0.4, 2)
  full <- lapply(d$rr.list, function(m) m | TRUE)
  empty <- lapply(d$rr.list, function(m) m & FALSE)
  p1 <- c(0.3, 0.7)
  p2 <- c(0.4, 0.2)
  expect_equal(bssr_power(full, d$rr.id, d$x11, d$x12, d$n21, d$n22, p1, p2,
                          d$n11, d$n12), c(1, 1), tolerance = 1e-12)
  expect_equal(bssr_power(empty, d$rr.id, d$x11, d$x12, d$n21, d$n22, p1, p2,
                          d$n11, d$n12), c(0, 0))
})

test_that("bssr_power rejects a region whose dimension does not fit the design", {
  d <- make_design(0.4, 3)
  expect_error(bssr_power(d$rr.list, d$rr.id, d$x11, d$x12, d$n21 + 1L, d$n22, 0.5, 0.5,
                          d$n11, d$n12), 'dimension')
  expect_error(bssr_power(d$rr.list, d$rr.id, d$x11, d$x12, d$n21[1], d$n22, 0.5, 0.5,
                          d$n11, d$n12), 'one element')
})
