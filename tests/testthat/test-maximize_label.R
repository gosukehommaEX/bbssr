test_that("maximize_label describes how the largest type I error rate was located", {
  x <- structure(list(), maximize = 'certified', interval = c(0.1, 0.9))
  expect_equal(maximize_label(x), 'certified over theta in [0.1, 0.9]')
  expect_equal(maximize_label(structure(list(), maximize = 'refined')),
               'largest values on the grid refined')
  expect_equal(maximize_label(structure(list(), maximize = 'grid')),
               'largest value on the grid')
  expect_equal(maximize_label(list()), 'not recorded')
})

test_that("maximize_label reports the intervals of the two boundaries", {
  x <- structure(list(), maximize = 'certified',
                 interval = matrix(c(0.1, 0.15, 0.9, 0.85), 2,
                                   dimnames = list(c('lower', 'upper'), c('from', 'to'))))
  expect_identical(maximize_label(x), paste('certified over theta in [0.1, 0.9] on the',
                                            'lower boundary and [0.15, 0.85] on the upper',
                                            'boundary'))
})

test_that("maximize_label names the nuisance parameter of the object", {
  x <- structure(list(), maximize = 'certified', interval = c(0.2, 1), nuisance = 'psi')
  expect_identical(maximize_label(x), 'certified over psi in [0.2, 1]')
  attr(x, 'nuisance') <- NULL
  expect_identical(maximize_label(x), 'certified over theta in [0.2, 1]')
})
