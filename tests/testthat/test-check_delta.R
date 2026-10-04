test_that("check_delta accepts effects in the direction of the alternative", {
  expect_null(check_delta(0.2, 'RD', 'greater'))
  expect_null(check_delta(-0.2, 'RD', 'less'))
  expect_null(check_delta(-0.2, 'RD', 'two.sided'))
  expect_null(check_delta(1.5, 'RR', 'greater'))
  expect_null(check_delta(0.5, 'OR', 'less'))
})

test_that("check_delta rejects invalid effects", {
  expect_error(check_delta(0, 'RD', 'greater'), 'differ from 0')
  expect_error(check_delta(1, 'OR', 'two.sided'), 'differ from 1')
  expect_error(check_delta(1.2, 'RD', 'greater'), 'risk difference')
  expect_error(check_delta(-0.5, 'RR', 'greater'), 'positive')
  expect_error(check_delta(-0.2, 'RD', 'greater'), 'exceed 0')
  expect_error(check_delta(2, 'RR', 'less'), 'fall below 1')
  expect_error(check_delta(c(0.1, 0.2), 'RD', 'greater', 'Delta.T'), 'Delta.T')
})
