test_that("check_delta accepts effects in the direction of the alternative", {
  expect_null(check_delta(0.2, 'RD', 'greater', margin = 0))
  expect_null(check_delta(-0.2, 'RD', 'less', margin = 0))
  expect_null(check_delta(-0.2, 'RD', 'two.sided', margin = 0))
  expect_null(check_delta(1.5, 'RR', 'greater', margin = 0))
  expect_null(check_delta(0.5, 'OR', 'less', margin = 0))
})

test_that("check_delta rejects invalid effects", {
  expect_error(check_delta(0, 'RD', 'greater', margin = 0), 'differ from 0')
  expect_error(check_delta(1, 'OR', 'two.sided', margin = 0), 'differ from 1')
  expect_error(check_delta(1.2, 'RD', 'greater', margin = 0), 'risk difference')
  expect_error(check_delta(-0.5, 'RR', 'greater', margin = 0), 'positive')
  expect_error(check_delta(-0.2, 'RD', 'greater', margin = 0), 'exceed 0')
  expect_error(check_delta(2, 'RR', 'less', margin = 0), 'fall below 1')
  expect_error(check_delta(c(0.1, 0.2), 'RD', 'greater', 'Delta.T', margin = 0), 'Delta.T')
})

test_that("check_delta accepts a margin with the risk difference and one-sided tests", {
  expect_null(check_delta(0, 'RD', 'greater', margin = 0.1))
  expect_null(check_delta(-0.05, 'RD', 'greater', margin = 0.1))
  expect_null(check_delta(0.05, 'RD', 'less', margin = 0.1))
  expect_error(check_delta(-0.1, 'RD', 'greater', margin = 0.1), 'exceed -margin')
  expect_error(check_delta(0.1, 'RD', 'less', margin = 0.1), 'fall below margin')
  expect_error(check_delta(1.5, 'RR', 'greater', margin = 0.1), "effect = 'RD'")
  expect_error(check_delta(0, 'RD', 'two.sided', margin = 0.1), 'one-sided')
  expect_error(check_delta(0, 'RD', 'greater', margin = NA), 'margin')
})
