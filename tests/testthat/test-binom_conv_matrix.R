test_that("binom_conv_matrix converts binomial probabilities to the Bernstein basis", {
  N <- 7
  M <- binom_conv_matrix(N, 0.2, 0.65)
  expect_equal(dim(M), c(N + 1, N + 1))
  expect_equal(colSums(M), rep(1, N + 1), tolerance = 1e-14)
  # b(X; N, a (1 - t) + c t) = sum_j M[X, j] B_{j,N}(t), with B_{j,N}(t) = b(j; N, t)
  for (t in c(0, 0.3, 0.8, 1)) {
    expect_equal(as.vector(M %*% dbinom(0:N, N, t)),
                 dbinom(0:N, N, 0.2 * (1 - t) + 0.65 * t), tolerance = 1e-14)
  }
})

test_that("each column of binom_conv_matrix is the distribution of a sum of binomials", {
  N <- 6
  M <- binom_conv_matrix(N, 0.1, 0.7)
  for (j in 0:N) {
    joint <- outer(dbinom(0:(N - j), N - j, 0.1), dbinom(0:j, j, 0.7))
    total <- outer(0:(N - j), 0:j, '+')
    expect_equal(M[, j + 1], as.vector(tapply(joint, total, sum)), tolerance = 1e-14)
  }
})

test_that("binom_conv_matrix is the identity on the unit interval", {
  expect_equal(binom_conv_matrix(5, 0, 1), diag(6), tolerance = 0)
  expect_equal(binom_conv_matrix(0, 0.3, 0.6), matrix(1), tolerance = 0)
})
