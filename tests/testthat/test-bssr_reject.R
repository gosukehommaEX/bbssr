test_that("bssr_reject passes the design to bssr_power", {
  map <- bssr_map(0.3, 10, 10, NULL, c(4, 3), 1, 0.025, 0.8, 'Chisq', FALSE, 'greater',
                  'minlike', 100L, 0, 'RD', 'standard', 'Chisq', 0.025, 'group', NULL,
                  NULL)
  st <- bssr_setup(map)
  rr.list <- lapply(seq_along(st$N1), function(k) {
    get_rr(st$N1[k], st$N2[k], 0.025, 'Chisq', 'greater', 'minlike', 100, 0)
  })
  got <- bssr_reject(st, rr.list, c(0.5, 0.3), c(0.2, 0.3))
  want <- bssr_power_ref(rr.list, st$rr.id, st$x11, st$x12, st$n21, st$n22,
                         c(0.5, 0.3), c(0.2, 0.3), st$n11, st$n12)
  expect_equal(got, want, tolerance = 1e-12)
})
