test_that("bern_restrict gives the coefficients on a subinterval", {
  set.seed(6)
  coef <- stats::runif(8)
  value <- function(b, t) sum(b * stats::dbinom(seq_along(b) - 1, length(b) - 1, t))
  for (iv in list(c(0.2, 0.7), c(0, 0.4), c(0.6, 1), c(0, 1))) {
    r <- bern_restrict(coef, iv[1], iv[2])
    for (s in c(0, 0.3, 0.8, 1)) {
      expect_equal(value(r, s), value(coef, iv[1] + (iv[2] - iv[1]) * s),
                   tolerance = 1e-13)
    }
  }
  expect_equal(bern_restrict(coef, 0, 1), coef)
})

test_that("bern_restrict returns a constant on a single point", {
  coef <- c(0.1, 0.5, 0.2, 0.9)
  v <- sum(coef * stats::dbinom(0:3, 3, 0.35))
  expect_equal(bern_restrict(coef, 0.35, 0.35), rep(v, 4))
  expect_equal(bern_restrict(2, 0.1, 0.4), 2)
})
