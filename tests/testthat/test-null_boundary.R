test_that("null_boundary returns the response probabilities on the null boundary", {
  b <- null_boundary(c(0.3, 0.5), 2, 'greater', 0.15)
  expect_equal(b$p1 - b$p2, c(-0.15, -0.15))
  expect_equal((2 * b$p1 + b$p2) / 3, c(0.3, 0.5))
  b <- null_boundary(0.5, 1, 'less', 0.1)
  expect_equal(c(b$p1, b$p2), c(0.55, 0.45))
})

test_that("null_boundary flags pooled probabilities without a valid boundary point", {
  b <- null_boundary(c(0, 0.02, 0.5, 0.99), 1, 'greater', 0.1)
  expect_identical(b$ok, c(FALSE, FALSE, TRUE, FALSE))
  expect_true(all(b$p1 >= 0 & b$p1 <= 1 & b$p2 >= 0 & b$p2 <= 1))
})

test_that("null_boundary returns theta for both groups without a margin", {
  theta <- c(0, 0.37, 1)
  b <- null_boundary(theta, 3, 'two.sided', 0)
  expect_identical(b$p1, theta)
  expect_identical(b$p2, theta)
  expect_true(all(b$ok))
})
