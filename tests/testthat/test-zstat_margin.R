test_that("zstat_margin reproduces the examples of Blackwelder (1982)", {
  # 30 patients per group. 18 responders with the standard and 13 with the experimental
  # therapy give z = 1.31 for the conventional hypothesis (p. 347). 18 responders with
  # the experimental and 21 with the standard therapy give 0.816 in absolute value for a
  # difference of 0.2 (p. 350), where the article subtracts in the other order
  z <- zstat_margin(30, 30, 0, 'unpooled')
  expect_equal(round(z[19, 14], 2), 1.31)
  expect_equal(z[19, 14], 1.31005615894438, tolerance = 1e-10)
  z <- zstat_margin(30, 30, 0.2, 'unpooled')
  expect_equal(round(z[19, 22], 3), 0.816)
  expect_equal(z[19, 22], 0.816496580927726, tolerance = 1e-10)
})

test_that("zstat_margin uses the restricted estimates for Farrington and Manning", {
  # Reference value from tools/reference/reference_values.py for 9 of 25 against 14 of 30
  # responders and a margin of 0.1
  z <- zstat_margin(25, 30, 0.1, 'restricted')
  expect_equal(z[10, 15], -0.050332228904682, tolerance = 1e-8)
  # Without a margin the restricted estimates are the pooled proportion
  expect_identical(zstat_margin(12, 9, 0, 'restricted'), zstat(12, 9))
})

test_that("zstat_margin handles a zero standard error", {
  z <- zstat_margin(10, 8, 0.2, 'unpooled')
  expect_identical(z[1, 1], Inf)
  expect_identical(z[11, 9], Inf)
  expect_identical(z[1, 9], -Inf)
  z <- zstat_margin(10, 8, 0, 'unpooled')
  expect_identical(z[1, 1], 0)
  expect_identical(z[11, 1], Inf)
  expect_identical(z[1, 9], -Inf)
  expect_true(all(is.finite(zstat_margin(10, 8, 0.2, 'restricted'))))
})
