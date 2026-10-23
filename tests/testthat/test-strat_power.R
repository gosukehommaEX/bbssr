test_that("strat_power agrees with the computation through the stored critical values", {
  # With two strata strat_fixed_power uses the stored critical values; the enumeration of
  # strat_power computes them as it goes
  p1 <- rbind(c(0.6, 0.45), c(0.2, 0.9))
  p2 <- rbind(c(0.35, 0.3), c(0.2, 0.9))
  for (type in 0:4) {
    for (alt in 0:2) {
      ts <- if (type == 4 && alt == 2) 1L else 0L
      a <- strat_power(c(8L, 6L), c(7L, 9L), p1, p2, type, alt, ts, 0.05)
      crit <- strat_crit(c(8L, 6L), c(7L, 9L), type, alt, ts, 0.05)
      b <- strat_bssr_power(c(0L, 0L), c(0L, 0L),
                            list(c(list(N1 = c(8L, 6L), N2 = c(7L, 9L)), crit)), 0L, p1, p2)
      expect_equal(a, b, tolerance = 1e-13, info = paste(type, alt))
    }
  }
})

test_that("strat_power checks its arguments", {
  expect_error(strat_power(c(4L, 5L), c(4L, 5L), matrix(0.5, 1, 3), matrix(0.5, 1, 3), 0L,
                           0L, 0L, 0.05), 'one column per stratum')
  expect_error(strat_power(c(4L, 5L), c(4L, 5L), matrix(1.5, 1, 2), matrix(0.5, 1, 2), 0L,
                           0L, 0L, 0.05), '\\[0, 1\\]')
})
