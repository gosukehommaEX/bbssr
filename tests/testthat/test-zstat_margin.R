test_that("zstat_margin reproduces the examples of Blackwelder (1982)", {
  # 30 patients per group. 18 responders with the standard and 13 with the experimental
  # therapy give z = 1.31 for the conventional hypothesis (p. 347). 18 responders with
  # the experimental and 21 with the standard therapy give 0.816 in absolute value for a
  # difference of 0.2 (p. 350), where the article subtracts in the other order
  z <- zstat_margin(30, 30, 0, 'unpooled', 'RD')
  expect_equal(round(z[19, 14], 2), 1.31)
  expect_equal(z[19, 14], 1.31005615894438, tolerance = 1e-10)
  z <- zstat_margin(30, 30, 0.2, 'unpooled', 'RD')
  expect_equal(round(z[19, 22], 3), 0.816)
  expect_equal(z[19, 22], 0.816496580927726, tolerance = 1e-10)
})

test_that("zstat_margin uses the restricted estimates for Farrington and Manning", {
  # Reference value from tools/reference/reference_values.py for 9 of 25 against 14 of 30
  # responders and a margin of 0.1
  z <- zstat_margin(25, 30, 0.1, 'restricted', 'RD')
  expect_equal(z[10, 15], -0.050332228904682, tolerance = 1e-8)
  # Without a margin the restricted estimates are the pooled proportion
  expect_identical(zstat_margin(12, 9, 0, 'restricted', 'RD'), zstat(12, 9))
})

test_that("zstat_margin handles a zero standard error", {
  z <- zstat_margin(10, 8, 0.2, 'unpooled', 'RD')
  expect_identical(z[1, 1], Inf)
  expect_identical(z[11, 9], Inf)
  expect_identical(z[1, 9], -Inf)
  z <- zstat_margin(10, 8, 0, 'unpooled', 'RD')
  expect_identical(z[1, 1], 0)
  expect_identical(z[11, 1], Inf)
  expect_identical(z[1, 9], -Inf)
  expect_true(all(is.finite(zstat_margin(10, 8, 0.2, 'restricted', 'RD'))))
})

test_that("zstat_margin computes the ratio statistic of Farrington and Manning (1990)", {
  # Reference values from tools/reference/reference_values.py for 15 of 25 against 20 of
  # 30 responders and the ratio 0.8, with the restricted and the observed proportions
  z.r <- zstat_margin(25, 30, 0.8, 'restricted', 'RR')
  z.u <- zstat_margin(25, 30, 0.8, 'unpooled', 'RR')
  expect_equal(c(z.r[16, 21], z.u[16, 21]), c(0.555141102984961, 0.556702214268904),
               tolerance = 1e-9)
  # The unpooled statistic written out directly
  se <- sqrt(0.6 * 0.4 / 25 + 0.8^2 * (20 / 30) * (10 / 30) / 30)
  expect_equal(z.u[16, 21], (0.6 - 0.8 * 20 / 30) / se, tolerance = 1e-12)
  # With the ratio 1 the restricted estimates are the pooled proportion
  expect_equal(zstat_margin(12, 9, 1, 'restricted', 'RR'), zstat(12, 9), tolerance = 1e-12)
})

test_that("zstat_margin handles a zero standard error of the ratio statistic", {
  z <- zstat_margin(10, 8, 0.8, 'unpooled', 'RR')
  expect_identical(z[1, 1], 0)
  expect_identical(z[11, 9], Inf)
  expect_identical(z[11, 1], Inf)
  expect_identical(z[1, 9], -Inf)
  z <- zstat_margin(10, 8, 0.8, 'restricted', 'RR')
  expect_true(all(is.finite(z)))
  expect_identical(z[1, 1], 0)
})
