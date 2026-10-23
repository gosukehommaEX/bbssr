test_that("paired_cp agrees with enumeration of the interim directions", {
  set.seed(4)
  N <- 12
  n1 <- 5
  R <- random_region(N)
  cp_ref <- function(D1, psi, pi) {
    n2 <- N - n1
    T2 <- trinomial_ref(n2, psi * pi, psi * (1 - pi))
    s <- 0
    for (a1 in 0:D1) {
      sub <- R[a1 + 1:(n2 + 1), D1 - a1 + 1:(n2 + 1)]
      s <- s + stats::dbinom(a1, D1, pi) * sum(T2 * sub)
    }
    s
  }
  D1 <- c(0L, 2L, 4L, 4L, 5L)
  psi <- c(0.3, 0.3, 0.3, 0.6, 0)
  pi <- c(0.7, 0.7, 0.7, 0.2, 0.5)
  expect_equal(paired_cp(R, n1, D1, psi, pi), mapply(cp_ref, D1, psi, pi),
               tolerance = 1e-12)
})

test_that("paired_cp reproduces the re-estimation of Schepers et al. (2026, Sec. 3.4)", {
  # 25 interim pairs with 5 discordant, the exact conditional McNemar test at the
  # two-sided level 0.05, psi = 0.2113 and p* = 0.8549. The article re-estimates 51
  # further pairs, so 76 pairs in all, and reports a conditional power of 0.78 with
  # psi = 0.2
  cp <- function(N, psi) {
    R <- get_paired_rr(N, 0.05, 'McNemar-exact', 'two.sided', 100, 0, FALSE, 0)
    paired_cp(R, 25L, 5L, psi, 0.854916)
  }
  expect_lt(cp(75, 0.211318), 0.8)
  expect_gte(cp(76, 0.211318), 0.8)
  expect_equal(round(cp(76, 0.2), 2), 0.78)
})

test_that("paired_cp checks its arguments", {
  R <- random_region(6)
  expect_error(paired_cp(R, 7L, 1L, 0.2, 0.5), 'n1 must lie')
  expect_error(paired_cp(R, 3L, 4L, 0.2, 0.5), 'D1 must lie')
  expect_error(paired_cp(R, 3L, c(1L, 2L), 0.2, 0.5), 'same length')
  expect_error(paired_cp(R, 3L, 1L, 1.2, 0.5), 'psi and pi_dir')
})
