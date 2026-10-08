test_that("null_boundary returns the response probabilities on the null boundary", {
  b <- null_boundary(c(0.3, 0.5), 2, 'greater', 0.15, 'RD')
  expect_equal(b$p1 - b$p2, c(-0.15, -0.15))
  expect_equal((2 * b$p1 + b$p2) / 3, c(0.3, 0.5))
  b <- null_boundary(0.5, 1, 'less', 0.1, 'RD')
  expect_equal(c(b$p1, b$p2), c(0.55, 0.45))
})

test_that("null_boundary flags pooled probabilities without a valid boundary point", {
  b <- null_boundary(c(0, 0.02, 0.5, 0.99), 1, 'greater', 0.1, 'RD')
  expect_identical(b$ok, c(FALSE, FALSE, TRUE, FALSE))
  expect_true(all(b$p1 >= 0 & b$p1 <= 1 & b$p2 >= 0 & b$p2 <= 1))
})

test_that("null_boundary returns theta for both groups without a margin", {
  theta <- c(0, 0.37, 1)
  b <- null_boundary(theta, 3, 'two.sided', 0, 'RD')
  expect_identical(b$p1, theta)
  expect_identical(b$p2, theta)
  expect_true(all(b$ok))
})

test_that("null_boundary follows the boundary of a ratio margin for both alternatives", {
  # With r = 2 and the ratio 0.8, p2 = 3 theta / 2.6 and p1 = 0.8 p2
  for (alt in c('greater', 'less')) {
    b <- null_boundary(c(0.3, 0.5, 0.9), 2, alt, 0.8, 'RR')
    expect_equal(b$p2[1:2], c(0.9, 1.5) / 2.6, info = alt)
    expect_equal(b$p1[1:2], 0.8 * b$p2[1:2], info = alt)
    expect_equal((2 * b$p1[1:2] + b$p2[1:2]) / 3, c(0.3, 0.5), info = alt)
    # At theta = 0.9, p2 = 2.7 / 2.6 exceeds 1 and is moved onto the unit interval, and
    # p1 = 0.8 x 2.7 / 2.6 is kept
    expect_identical(b$ok, c(TRUE, TRUE, FALSE), info = alt)
    expect_equal(c(b$p1[3], b$p2[3]), c(0.8 * 2.7 / 2.6, 1), info = alt)
  }
})
