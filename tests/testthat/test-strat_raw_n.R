test_that("strat_raw_n agrees with the reference values of Woolson et al. and Nam", {
  expect_equal(c(strat_raw_n(c(0.5, 0.7), c(0.3, 0.5), c(0.4, 0.6), c(1, 1), 0.025, 0.8,
                             'greater', FALSE),
                 strat_raw_n(c(0.5, 0.7), c(0.3, 0.5), c(0.4, 0.6), c(1, 1), 0.025, 0.8,
                             'greater', TRUE),
                 strat_raw_n(c(0.3, 0.5, 0.2), c(0.45, 0.6, 0.35), c(1, 2, 1) / 4,
                             c(2, 2, 2), 0.05, 0.9, 'less', FALSE),
                 strat_raw_n(c(0.6, 0.4), c(0.4, 0.3), c(0.5, 0.5), c(0.5, 2), 0.05, 0.8,
                             'two.sided', TRUE)),
               c(185.997689655091, 205.511097928377, 566.558428740696, 402.041885184174),
               tolerance = 1e-10)
})

test_that("strat_raw_n is vectorized over scenarios and NaN without an effect", {
  p1 <- rbind(c(0.5, 0.7), c(0.3, 0.5), c(0.6, 0.4))
  p2 <- rbind(c(0.3, 0.5), c(0.3, 0.5), c(0.4, 0.6))
  n <- strat_raw_n(p1, p2, c(0.5, 0.5), c(1, 1), 0.025, 0.8, 'greater', FALSE)
  expect_length(n, 3)
  expect_true(is.finite(n[1]))
  # No effect, and effects that cancel over the strata
  expect_true(is.nan(n[2]) && is.nan(n[3]))
  expect_true(is.nan(strat_raw_n(c(0.5, 0.7), c(0.3, 0.5), c(0.5, 0.5), c(1, 1), 0.025,
                                 0.8, 'less', FALSE)))
  # With one stratum the formula is that of two independent binomial proportions
  one <- strat_raw_n(c(0.5, 0.5), c(0.3, 0.3), c(0.5, 0.5), c(1, 1), 0.025, 0.8,
                     'greater', FALSE)
  z <- stats::qnorm(c(0.975, 0.8))
  ref <- 2 * (z[1] * sqrt(2 * 0.4 * 0.6) + z[2] * sqrt(0.5 * 0.5 + 0.3 * 0.7))^2 / 0.2^2
  expect_equal(one, ref)
})
