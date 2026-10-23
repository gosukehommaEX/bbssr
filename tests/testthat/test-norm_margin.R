test_that("norm_margin returns the margins of a test of equivalence", {
  expect_identical(norm_margin('greater', 0.1, 'RD'), 0.1)
  expect_identical(norm_margin('two.sided', 0, 'RD'), 0)
  expect_equal(norm_margin('equivalence', 0.15, 'RD'), c(0.15, 0.15))
  expect_equal(norm_margin('equivalence', c(0.1, 0.2), 'RD'), c(0.1, 0.2))
  expect_equal(norm_margin('equivalence', 1.25, 'RR'), c(0.8, 1.25))
  expect_equal(norm_margin('equivalence', 0.8, 'RR'), c(0.8, 1.25))
  expect_equal(norm_margin('equivalence', c(0.75, 1.4), 'RR'), c(0.75, 1.4))
})

test_that("norm_margin rejects invalid margins of a test of equivalence", {
  expect_error(norm_margin('equivalence', c(0.1, 0.2, 0.3), 'RD'), 'one or two')
  expect_error(norm_margin('equivalence', NA_real_, 'RD'), 'one or two')
  expect_error(norm_margin('equivalence', 0, 'RD'), 'in \\(0, 1\\)')
  expect_error(norm_margin('equivalence', c(-0.1, 0.2), 'RD'), 'in \\(0, 1\\)')
  expect_error(norm_margin('equivalence', 1, 'RR'), 'differ from 1')
  expect_error(norm_margin('equivalence', c(1.1, 1.25), 'RR'), '0 < lower < 1 < upper')
})
