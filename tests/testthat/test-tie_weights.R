# Design of test-bssr_reject.R: Delta.A = 0.3, N1 = N2 = 10, interim sizes 4 and 3
weights_design <- function() {
  map <- bssr_map(0.3, 10, 10, NULL, c(4, 3), 1, 0.025, 0.8, 'Chisq', FALSE, 'greater',
                  'minlike', 100L, 0, 'RD', 'standard', 'Chisq', 0.025, 'group', NULL,
                  NULL, FALSE, 0)
  st <- bssr_setup(map)
  rr.list <- lapply(seq_along(st$N1), function(k) {
    get_rr(st$N1[k], st$N2[k], 0.025, 'Chisq', 'greater', 'minlike', 100, 0, FALSE, 0)
  })
  list(st = st, rr.list = rr.list)
}

# Rejection probability from the weights of the final responder counts
weights_value <- function(w, p1, p2) {
  sum(vapply(w, function(x) {
    sum(dbinom(0:x$N1, x$N1, p1) * (x$C %*% dbinom(0:x$N2, x$N2, p2)))
  }, numeric(1)))
}

test_that("tie_weights reproduces the rejection probability of a re-estimation design", {
  d <- weights_design()
  w <- tie_weights(d$st, d$rr.list)
  expect_length(w, length(d$st$N1))
  p1 <- c(0.5, 0.3, 0.12)
  p2 <- c(0.2, 0.3, 0.4)
  got <- vapply(1:3, function(i) weights_value(w, p1[i], p2[i]), numeric(1))
  expect_equal(got, bssr_reject(d$st, d$rr.list, p1, p2), tolerance = 1e-13)
})

test_that("the weights vanish outside the rejection region and add up to one", {
  d <- weights_design()
  w <- tie_weights(d$st, d$rr.list)
  for (k in seq_along(w)) {
    expect_equal(dim(w[[k]]$C), c(d$st$N1[k] + 1, d$st$N2[k] + 1))
    expect_true(all(w[[k]]$C[!d$rr.list[[k]]] == 0))
  }
  # With every outcome rejected the rejection probability is one
  everything <- lapply(d$rr.list, function(m) matrix(TRUE, nrow(m), ncol(m)))
  expect_equal(weights_value(tie_weights(d$st, everything), 0.37, 0.61), 1, tolerance = 1e-14)
})

test_that("tie_weights of a fixed-sample design is its rejection region", {
  rr <- get_rr(12, 9, 0.025, 'Chisq', 'greater', 'minlike', 100, 0, FALSE, 0)
  w <- tie_weights(fixed_setup(12, 9), list(rr))
  expect_equal(w[[1]]$C, rr * 1, tolerance = 0)
})
