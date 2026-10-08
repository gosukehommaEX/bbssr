test_that("format_margin labels the margin of a result for printing", {
  expect_null(format_margin(NULL, NULL))
  expect_null(format_margin(0, NULL))
  expect_null(format_margin(0, 'RD'))
  expect_identical(format_margin(0.1, NULL), '0.1')
  expect_identical(format_margin(-0.05, 'RD'), '-0.05')
  expect_identical(format_margin(0.8, 'RR'), '0.8 (risk ratio)')
  expect_identical(format_margin(1, 'RR'), '1 (risk ratio)')
})
