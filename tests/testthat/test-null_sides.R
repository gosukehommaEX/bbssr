test_that("null_sides keeps a one-sided or two-sided null hypothesis whole", {
  s <- null_sides('greater', 0.1, 'RD')
  expect_length(s, 1)
  expect_identical(s[[1]]$alternative, 'greater')
  expect_identical(s[[1]]$margin, 0.1)
  expect_true(is.na(s[[1]]$side))
  expect_identical(null_sides('two.sided', 0, 'RD')[[1]]$alternative, 'two.sided')
})

test_that("null_sides splits a test of equivalence into two one-sided components", {
  s <- null_sides('equivalence', c(0.1, 0.2), 'RD')
  expect_length(s, 2)
  expect_identical(vapply(s, `[[`, '', 'alternative'), c('greater', 'less'))
  expect_identical(vapply(s, `[[`, 0, 'margin'), c(0.1, 0.2))
  expect_identical(vapply(s, `[[`, '', 'side'), c('lower', 'upper'))
  s <- null_sides('equivalence', 1.25, 'RR')
  expect_equal(vapply(s, `[[`, 0, 'margin'), c(0.8, 1.25))
})
