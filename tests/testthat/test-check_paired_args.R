test_that("check_paired_args returns the arguments in computational form", {
  a <- check_paired_args(20, 0.05, 'McNemar-u', 50, 0.001, TRUE, 'greater', 0)
  expect_identical(a$N, 20L)
  expect_identical(a$n.grid, 50L)
  expect_identical(a$Test, 'McNemar-uncond')
  expect_true(a$ref.pvalue)
  # The refinement is switched off for the tests without a nuisance parameter
  expect_false(check_paired_args(20, 0.05, 'McNemar-exact', 50, 0, TRUE, 'greater',
                                 0)$ref.pvalue)
  # A single margin of equivalence stands for both
  expect_equal(check_paired_args(20, 0.05, 'Paired-Score', 50, 0, FALSE, 'equivalence',
                                 0.1)$margin, c(0.1, 0.1))
})

test_that("check_paired_args rejects invalid arguments", {
  ok <- function(...) {
    args <- list(N = 20, alpha = 0.05, Test = 'McNemar', n.grid = 50, bb.gamma = 0,
                 ref.pvalue = FALSE, alternative = 'greater', margin = 0)
    new <- list(...)
    args[names(new)] <- new
    do.call(check_paired_args, args)
  }
  expect_error(ok(Test = 'Chisq'), "should be one of")
  expect_error(ok(N = 2.5), 'positive integer')
  expect_error(ok(alpha = 1), 'alpha')
  expect_error(ok(n.grid = 1), 'n.grid')
  expect_error(ok(bb.gamma = 0.06), 'bb.gamma')
  expect_error(ok(ref.pvalue = NA), 'ref.pvalue')
  expect_error(ok(margin = 0.1), "requires Test = 'Paired-Score'")
  expect_error(ok(Test = 'Paired-Score', margin = 0.1, alternative = 'two.sided'),
               'one-sided')
  expect_error(ok(margin = 1.2), 'margin must be')
  expect_error(ok(alternative = 'equivalence', margin = 0.1),
               "'equivalence' requires Test = 'Paired-Score'")
  expect_error(ok(Test = 'Paired-Score', alternative = 'equivalence', margin = -0.1),
               'in \\(0, 1\\)')
  expect_warning(ok(bb.gamma = 0.001), 'only to the McNemar-uncond')
})
