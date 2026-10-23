test_that("paired_bssr_power agrees with enumeration of both stages", {
  set.seed(2)
  n1 <- 5
  sizes <- c(8, 8, 10, 10, 12, 12)
  Ns <- c(8, 10, 12)
  regions <- lapply(Ns, random_region)
  names(regions) <- Ns
  rr.id <- match(sizes, Ns) - 1L
  for (p in list(c(0.3, 0.1), c(0.25, 0.25), c(0, 0.4), c(0.6, 0.4), c(0, 0))) {
    expect_equal(paired_bssr_power(unname(regions), rr.id, n1, p[1], p[2]),
                 paired_bssr_ref(sizes, regions, n1, p[1], p[2]), tolerance = 1e-12,
                 info = paste(p, collapse = ', '))
  }
})

test_that("paired_bssr_power reduces to paired_power without re-estimation", {
  set.seed(3)
  R <- random_region(10)
  expect_equal(paired_bssr_power(list(R), rep(0L, 5), 4, c(0.2, 0.4), c(0.1, 0.3)),
               paired_power(R, c(0.2, 0.4), c(0.1, 0.3)), tolerance = 1e-13)
})

test_that("paired_bssr_power checks its arguments", {
  R <- random_region(6)
  expect_error(paired_bssr_power(list(R), rep(0L, 3), 3, 0.1, 0.1), 'n1 \\+ 1')
  expect_error(paired_bssr_power(list(R), rep(1L, 4), 3, 0.1, 0.1), 'positions')
  expect_error(paired_bssr_power(list(R), rep(0L, 9), 8, 0.1, 0.1), 'below the interim')
  expect_error(paired_bssr_power(list(R + 0L), rep(0L, 4), 3, 0.1, 0.1),
               'logical matrices')
})
