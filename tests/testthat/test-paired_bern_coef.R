test_that("paired_bern_coef gives the rejection probability on each boundary", {
  set.seed(5)
  N <- 9
  R <- random_region(N)
  W <- stats::runif(N + 1)
  D <- outer(0:N, 0:N, '+')
  direct <- function(p10, p01) {
    P <- trinomial_ref(N, p10, p01)
    sum(P * R * W[pmin(D, N) + 1])
  }
  for (case in list(c(0, 0), c(1, 0.15), c(2, 0.15), c(1, 0.4))) {
    beta <- paired_bern_coef(R, W, case[1], case[2])
    m <- case[2]
    for (t in c(0, 0.17, 0.5, 0.83, 1)) {
      psi <- if (case[1] == 0) t else m + (1 - m) * t
      d0 <- c(0, -m, m)[case[1] + 1]
      value <- sum(beta * stats::dbinom(0:N, N, t))
      expect_equal(value, direct((psi + d0) / 2, (psi - d0) / 2), tolerance = 1e-12,
                   info = paste(c(case, t), collapse = ', '))
    }
  }
})

test_that("paired_bern_coef checks its arguments", {
  R <- random_region(4)
  expect_error(paired_bern_coef(R, rep(1, 4), 0, 0), 'N \\+ 1')
  expect_error(paired_bern_coef(R, rep(1, 5), 3, 0), 'type')
  expect_error(paired_bern_coef(R, rep(1, 5), 1, 1), 'm must lie')
})
