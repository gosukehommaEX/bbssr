test_that("null_range spans theta without a margin", {
  expect_equal(null_range(c(0.3, 0.1, 0.7), 1, 'greater', 0, 'RD'), c(0.1, 0.7))
  expect_equal(null_range(seq(0, 1, by = 0.005), 2, 'two.sided', 0, 'RD'), c(0, 1))
})

test_that("null_range stops where a response probability leaves the unit interval", {
  theta <- seq(0, 1, by = 0.005)
  # With r = 2 and a margin of 0.2, p1 = theta - 0.2 / 3 and p2 = theta + 0.4 / 3
  iv <- null_range(theta, 2, 'greater', 0.2, 'RD')
  expect_equal(iv, c(0.2 / 3, 1 - 0.4 / 3), tolerance = 1e-14)
  b <- null_boundary(iv, 2, 'greater', 0.2, 'RD')
  expect_equal(c(b$p1[1], b$p2[2]), c(0, 1), tolerance = 1e-14)
  # For alternative = 'less' the boundary is p1 - p2 = 0.2
  expect_equal(null_range(theta, 2, 'less', 0.2, 'RD'), c(0.4 / 3, 1 - 0.2 / 3),
               tolerance = 1e-14)
  # The interval does not extend beyond the values of theta given
  expect_equal(null_range(c(0.05, 0.3, 0.5, 0.95), 1, 'greater', 0.2, 'RD'), c(0.1, 0.9))
  expect_equal(null_range(c(0.3, 0.5, 0.7), 1, 'greater', 0.2, 'RD'), c(0.3, 0.7))
})

test_that("null_range keeps a value of theta accepted within rounding error", {
  theta <- c(0.1 - 1e-12, 0.5)
  expect_true(all(null_boundary(theta, 1, 'greater', 0.2, 'RD')$ok))
  expect_identical(null_range(theta, 1, 'greater', 0.2, 'RD'), theta)
})

test_that("null_range stops where a response probability on a ratio boundary reaches 1", {
  theta <- seq(0, 1, by = 0.005)
  # With r = 2 and the ratio 0.8, p2 = 3 theta / 2.6 reaches 1 at theta = 2.6 / 3
  iv <- null_range(theta, 2, 'greater', 0.8, 'RR')
  expect_equal(iv, c(0, 2.6 / 3), tolerance = 1e-14)
  b <- null_boundary(iv, 2, 'greater', 0.8, 'RR')
  expect_equal(c(b$p1[1], b$p2[1], b$p2[2]), c(0, 0, 1), tolerance = 1e-14)
  # With the ratio 1.5, p1 = 1.5 p2 reaches 1 first, at theta = 4 / 4.5
  iv <- null_range(theta, 2, 'less', 1.5, 'RR')
  expect_equal(iv, c(0, 4 / 4.5), tolerance = 1e-14)
  expect_equal(null_boundary(iv[2], 2, 'less', 1.5, 'RR')$p1, 1, tolerance = 1e-14)
  expect_equal(null_range(c(0.2, 0.6), 1, 'greater', 0.5, 'RR'), c(0.2, 0.6))
})
