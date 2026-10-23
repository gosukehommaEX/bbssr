test_that("ss_raw_n2_equiv matches independent reference values", {
  # Reference values from tools/reference/reference_values.py, whose restricted estimates
  # are found by bisection on the score
  n2 <- c(ss_raw_n2_equiv(0.5, 0.5, 1, 0.05, 0.8, 'standard', c(0.15, 0.15), 'RD'),
          ss_raw_n2_equiv(0.3, 0.35, 2, 0.025, 0.9, 'alternative.variance', c(0.1, 0.2),
                          'RD'),
          ss_raw_n2_equiv(0.3, 0.35, 2, 0.025, 0.9, 'null.variance', c(0.1, 0.2), 'RD'),
          ss_raw_n2_equiv(0.6, 0.6, 1, 0.05, 0.8, 'standard', c(0.8, 1.25), 'RR'))
  expect_equal(n2, c(187.894968505217, 1397.4872671716, 1423.30926399004,
                     234.224190102778), tolerance = 1e-8)
  # ss_raw_n2 hands a test of equivalence to it
  expect_equal(ss_raw_n2(0.5, 0.5, 1, 0.05, 0.8, 'equivalence', 'standard', 0.15, 'RD'),
               n2[1], tolerance = 1e-12)
})

test_that("the approximate power of ss_raw_n2_equiv reaches the target at the root", {
  lo <- ss_raw_parts(0.4, 0.45, 1, 'greater', 'standard', 0.2, 'RD')
  up <- ss_raw_parts(0.4, 0.45, 1, 'less', 'standard', 0.2, 'RD')
  n <- ss_raw_n2_equiv(0.4, 0.45, 1, 0.05, 0.85, 'standard', 0.2, 'RD')
  za <- qnorm(0.95)
  pw <- pnorm((lo$d * sqrt(n) - za * sqrt(lo$v0)) / sqrt(lo$v1)) +
    pnorm((up$d * sqrt(n) - za * sqrt(up$v0)) / sqrt(up$v1)) - 1
  expect_equal(pw, 0.85, tolerance = 1e-9)
  # Vectorized over the probabilities, NaN outside the equivalence region
  v <- ss_raw_n2_equiv(c(0.4, 0.4, 0.4), c(0.45, 0.65, 0.4), 1, 0.05, 0.85, 'standard',
                       0.2, 'RD')
  expect_equal(v[1], n, tolerance = 1e-12)
  expect_true(is.nan(v[2]))
  expect_true(is.finite(v[3]))
  expect_true(is.na(ss_raw_n2_equiv(NA_real_, 0.4, 1, 0.05, 0.85, 'standard', 0.2, 'RD')))
  # With alpha above 0.5 the approximation exceeds the target without patients
  expect_identical(ss_raw_n2_equiv(0.4, 0.4, 1, 0.9, 0.2, 'standard', 0.2, 'RD'), 0)
})
