test_that("check_delta accepts effects in the direction of the alternative", {
  expect_null(check_delta(0.2, 'RD', 'greater', margin = 0, margin.scale = 'RD'))
  expect_null(check_delta(-0.2, 'RD', 'less', margin = 0, margin.scale = 'RD'))
  expect_null(check_delta(-0.2, 'RD', 'two.sided', margin = 0, margin.scale = 'RD'))
  expect_null(check_delta(1.5, 'RR', 'greater', margin = 0, margin.scale = 'RD'))
  expect_null(check_delta(0.5, 'OR', 'less', margin = 0, margin.scale = 'RD'))
})

test_that("check_delta rejects invalid effects", {
  expect_error(check_delta(0, 'RD', 'greater', margin = 0, margin.scale = 'RD'),
               'differ from 0')
  expect_error(check_delta(1, 'OR', 'two.sided', margin = 0, margin.scale = 'RD'),
               'differ from 1')
  expect_error(check_delta(1.2, 'RD', 'greater', margin = 0, margin.scale = 'RD'),
               'risk difference')
  expect_error(check_delta(-0.5, 'RR', 'greater', margin = 0, margin.scale = 'RD'),
               'positive')
  expect_error(check_delta(-0.2, 'RD', 'greater', margin = 0, margin.scale = 'RD'),
               'exceed 0')
  expect_error(check_delta(2, 'RR', 'less', margin = 0, margin.scale = 'RD'),
               'fall below 1')
  expect_error(check_delta(c(0.1, 0.2), 'RD', 'greater', 'Delta.T', margin = 0,
                           margin.scale = 'RD'), 'Delta.T')
})

test_that("check_delta accepts a margin with the risk difference and one-sided tests", {
  expect_null(check_delta(0, 'RD', 'greater', margin = 0.1, margin.scale = 'RD'))
  expect_null(check_delta(-0.05, 'RD', 'greater', margin = 0.1, margin.scale = 'RD'))
  expect_null(check_delta(0.05, 'RD', 'less', margin = 0.1, margin.scale = 'RD'))
  expect_error(check_delta(-0.1, 'RD', 'greater', margin = 0.1, margin.scale = 'RD'),
               'exceed -margin')
  expect_error(check_delta(0.1, 'RD', 'less', margin = 0.1, margin.scale = 'RD'),
               'fall below margin')
  expect_error(check_delta(1.5, 'RR', 'greater', margin = 0.1, margin.scale = 'RD'),
               "effect = 'RD'")
  expect_error(check_delta(0, 'RD', 'two.sided', margin = 0.1, margin.scale = 'RD'),
               'one-sided')
  expect_error(check_delta(0, 'RD', 'greater', margin = NA, margin.scale = 'RD'),
               'margin')
})

test_that("check_delta accepts a ratio margin with the risk ratio and one-sided tests", {
  expect_null(check_delta(1, 'RR', 'greater', margin = 0.8, margin.scale = 'RR'))
  expect_null(check_delta(0.9, 'RR', 'greater', margin = 0.8, margin.scale = 'RR'))
  expect_null(check_delta(1, 'RR', 'less', margin = 1.5, margin.scale = 'RR'))
  expect_error(check_delta(0.8, 'RR', 'greater', margin = 0.8, margin.scale = 'RR'),
               'exceed margin')
  expect_error(check_delta(1.5, 'RR', 'less', margin = 1.5, margin.scale = 'RR'),
               'fall below margin')
  expect_error(check_delta(0, 'RD', 'greater', margin = 0.8, margin.scale = 'RR'),
               "effect = 'RR'")
  expect_error(check_delta(1, 'RR', 'two.sided', margin = 0.8, margin.scale = 'RR'),
               'one-sided')
  expect_error(check_delta(-1, 'RR', 'less', margin = 0.8, margin.scale = 'RR'),
               'positive')
  for (m in list(0, -0.5, NA, Inf, c(0.8, 0.9))) {
    expect_error(check_delta(1, 'RR', 'greater', margin = m, margin.scale = 'RR'),
                 'positive value', info = paste(m, collapse = ' '))
  }
})
