test_that("check_strat_args resolves the test and returns the codes of the compiled code", {
  a <- check_strat_args(c(5, 6), c(4, 7), 0.05, 'Strat', 'two.sided', 'blaker')
  expect_identical(a$Test, 'Strat-Fisher')
  expect_identical(c(a$type, a$alt, a$ts), c(4L, 2L, 2L))
  expect_identical(a$N1, c(5L, 6L))
  expect_false(a$cc)
  b <- check_strat_args(c(5, 6), c(4, 7), 0.05, 'Cochran-cc', 'less', 'blaker')
  # The two-sided method matters only for the stratified exact test with a two-sided
  # alternative, and is set to minlike otherwise
  expect_identical(c(b$type, b$alt, b$ts), c(3L, 1L, 0L))
  expect_identical(b$tsmethod, 'minlike')
  expect_true(b$cc)
  expect_identical(check_strat_args(c(5, 6), c(4, 7), 0.05, 'Cochran', 'greater',
                                    'minlike')$Test, 'Cochran')
})

test_that("check_strat_args rejects invalid designs", {
  expect_error(check_strat_args(5, 4, 0.05, 'CMH', 'greater', 'minlike'), 'at least')
  expect_error(check_strat_args(c(5, 6), c(4, 7, 3), 0.05, 'CMH', 'greater', 'minlike'),
               'same length')
  expect_error(check_strat_args(c(5, 0), c(4, 7), 0.05, 'CMH', 'greater', 'minlike'),
               'positive integers')
  expect_error(check_strat_args(c(5, 6.5), c(4, 7), 0.05, 'CMH', 'greater', 'minlike'),
               'positive integers')
  expect_error(check_strat_args(c(5, 6), c(4, 7), 1, 'CMH', 'greater', 'minlike'), 'alpha')
  expect_error(check_strat_args(c(5, 6), c(4, 7), 0.05, 'Chisq', 'greater', 'minlike'))
  expect_error(check_strat_args(c(5, 6), c(4, 7), 0.05, 'CMH', 'equivalence', 'minlike'))
})

test_that("check_strat_args returns the codes and flags of the unstratified tests", {
  u <- function(Test, alternative = 'two.sided') {
    check_strat_args(c(5, 6), c(4, 7), 0.05, Test, alternative, 'blaker')
  }
  expect_identical(c(u('Unstrat-Chisq')$type, u('Unstrat-Chisq-cc')$type,
                     u('Unstrat-Fisher')$type), c(5L, 6L, 7L))
  expect_identical(c(u('Unstrat-Chisq')$cc, u('Unstrat-Chisq-cc')$cc), c(FALSE, TRUE))
  expect_identical(c(u('Unstrat-Chisq')$unstrat, u('Unstrat-Fisher')$unstrat,
                     u('Strat-Fisher')$unstrat, u('CMH')$unstrat),
                   c(TRUE, TRUE, FALSE, FALSE))
  # The two-sided method is kept for the two exact tests only
  expect_identical(c(u('Unstrat-Fisher')$tsmethod, u('Unstrat-Chisq')$tsmethod),
                   c('blaker', 'minlike'))
  expect_identical(u('Unstrat-Fisher', 'greater')$tsmethod, 'minlike')
  expect_error(u('Unstrat'))
})

test_that("check_strat_args refuses infinite and too large sizes", {
  expect_error(check_strat_args(c(5, Inf), c(4, 7), 0.05, 'CMH', 'greater', 'minlike'),
               'positive integers')
  expect_error(check_strat_args(c(5, 6), c(4, 2^31), 0.05, 'CMH', 'greater', 'minlike'),
               'positive integers')
})
