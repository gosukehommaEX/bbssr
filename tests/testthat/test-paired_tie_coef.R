test_that("the Bernstein coefficients give the rejection probability of the design", {
  set.seed(7)
  n1 <- 4
  N.D1 <- c(8, 8, 10, 12, 12)
  Ns <- c(8, 10, 12)
  rr <- lapply(Ns, random_region)
  W <- paired_tie_weights(N.D1, n1, Ns)
  rr.id <- match(N.D1, Ns) - 1L
  for (case in list(c('greater', 0), c('greater', 0.2), c('less', 0.1))) {
    b0 <- paired_boundary(0.5, case[1], as.numeric(case[2]))
    coef <- paired_tie_coef(rr, Ns, W, b0$type, b0$m)
    expect_length(coef, 13)
    for (psi in c(0.25, 0.5, 0.9)) {
      b <- paired_boundary(psi, case[1], as.numeric(case[2]))
      expect_equal(sum(coef * stats::dbinom(0:12, 12, b$t)),
                   paired_bssr_power(rr, rr.id, n1, b$p10, b$p01), tolerance = 1e-12,
                   info = paste(c(case, psi), collapse = ', '))
    }
  }
})
