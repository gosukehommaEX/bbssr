test_that("fixed_setup describes a design without an interim stage", {
  st <- fixed_setup(12, 9)
  expect_identical(c(st$n11, st$n12, st$x11, st$x12, st$rr.id), rep(0L, 5))
  expect_equal(c(st$N1, st$N2, st$n21, st$n22), c(12, 9, 12, 9))
  # The rejection probability through bssr_reject is that of the fixed-sample design
  rr <- get_rr(12, 9, 0.025, 'Chisq', 'greater', 'minlike', 100, 0, FALSE, 0, 'RD')
  p1 <- c(0.2, 0.5, 0.7)
  p2 <- c(0.2, 0.3, 0.6)
  want <- vapply(1:3, function(i) {
    power_from_rr(rr, dbinom(0:12, 12, p1[i]), dbinom(0:9, 9, p2[i]))
  }, numeric(1))
  expect_equal(bssr_reject(st, list(rr), p1, p2), want, tolerance = 1e-13)
})
