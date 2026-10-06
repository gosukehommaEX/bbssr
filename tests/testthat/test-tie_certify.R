# Design of test-bssr_reject.R: Delta.A = 0.3, N1 = N2 = 10, interim sizes 4 and 3
certify_design <- function() {
  map <- bssr_map(0.3, 10, 10, NULL, c(4, 3), 1, 0.025, 0.8, 'Chisq', FALSE, 'greater',
                  'minlike', 100L, 0, 'RD', 'standard', 'Chisq', 0.025, 'group', NULL,
                  NULL, FALSE, 0)
  st <- bssr_setup(map)
  rr.list <- lapply(seq_along(st$N1), function(k) {
    get_rr(st$N1[k], st$N2[k], 0.025, 'Chisq', 'greater', 'minlike', 100, 0, FALSE, 0)
  })
  list(st = st, rr.list = rr.list)
}

test_that("tie_certify bounds the type I error rate over the interval", {
  d <- certify_design()
  w <- tie_weights(d$st, d$rr.list)
  for (iv in list(c(0, 1), c(0.2, 0.6))) {
    cm <- tie_certify(w, iv, 1, 'greater', 0)
    th <- seq(iv[1], iv[2], length.out = 2001)
    expect_gte(cm$y, max(bssr_reject(d$st, d$rr.list, th, th)) - 1e-12)
    expect_gte(cm$bound, cm$y)
    expect_lte(cm$bound - cm$y, 1e-12)
    expect_true(cm$x >= iv[1] && cm$x <= iv[2])
    expect_equal(bssr_reject(d$st, d$rr.list, cm$x, cm$x), cm$y, tolerance = 1e-12)
  }
})

test_that("tie_certify follows the boundary of a non-inferiority hypothesis", {
  rr <- get_rr(20, 20, 0.025, 'Farrington-Manning', 'greater', 'minlike', 100, 0, FALSE,
               0.2)
  w <- tie_weights(fixed_setup(20, 20), list(rr))
  iv <- null_range(seq(0, 1, by = 0.01), 1, 'greater', 0.2)
  cm <- tie_certify(w, iv, 1, 'greater', 0.2)
  rate <- function(theta) {
    b <- null_boundary(theta, 1, 'greater', 0.2)
    vapply(seq_along(theta), function(i) {
      power_from_rr(rr, dbinom(0:20, 20, b$p1[i]), dbinom(0:20, 20, b$p2[i]))
    }, numeric(1))
  }
  expect_gte(cm$y, max(rate(seq(iv[1], iv[2], length.out = 2001))) - 1e-12)
  expect_lte(cm$bound - cm$y, 1e-12)
  expect_equal(rate(cm$x), cm$y, tolerance = 1e-12)
})
