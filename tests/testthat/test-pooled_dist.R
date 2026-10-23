test_that("pooled_dist is the convolution of the two binomial distributions", {
  d <- pooled_dist(4, 6, 0.3, 0.55)
  expect_length(d, 11)
  expect_equal(sum(d), 1)
  ref <- vapply(0:10, function(s) {
    sum(stats::dbinom(0:4, 4, 0.3) * stats::dbinom(s - 0:4, 6, 0.55))
  }, numeric(1))
  expect_equal(d, ref, tolerance = 1e-15)
})
