test_that("check_rr_args returns the arguments in computational form", {
  a <- check_rr_args(10, 8, 0.025, 'Bosch', 50, 0, TRUE)
  expect_identical(a, list(N1 = 10L, N2 = 8L, n.grid = 50L, Test = 'Boschloo',
                           ref.pvalue = TRUE))
})

test_that("check_rr_args rejects invalid arguments", {
  expect_error(check_rr_args(c(5, 6), 5, 0.05, 'Chisq', 100, 0, FALSE), 'single value')
  expect_error(check_rr_args(0, 5, 0.05, 'Chisq', 100, 0, FALSE), 'positive integers')
  expect_error(check_rr_args(5, 4.5, 0.05, 'Chisq', 100, 0, FALSE), 'positive integers')
  expect_error(check_rr_args(5, 5, 0, 'Chisq', 100, 0, FALSE), 'alpha')
  expect_error(check_rr_args(5, 5, 0.05, 'Chisq', 1, 0, FALSE), 'n.grid')
  expect_error(check_rr_args(5, 5, 0.05, 'Z-pool', 100, 0.05, FALSE), 'bb.gamma')
  expect_error(check_rr_args(5, 5, 0.05, 'Student', 100, 0, FALSE), 'should be one of')
})

test_that("check_rr_args warns when bb.gamma is given for a conditional test", {
  expect_warning(check_rr_args(5, 5, 0.05, 'Fisher', 100, 0.001, FALSE), 'bb.gamma')
  expect_no_warning(check_rr_args(5, 5, 0.05, 'Z-pool', 100, 0.001, FALSE))
})

test_that("check_rr_args validates ref.pvalue and keeps it for unconditional tests", {
  expect_true(check_rr_args(5, 5, 0.05, 'Z-pool', 100, 0, TRUE)$ref.pvalue)
  expect_true(check_rr_args(5, 5, 0.05, 'Boschloo', 100, 0, TRUE)$ref.pvalue)
  for (tst in c('Chisq', 'Fisher', 'Fisher-midP')) {
    expect_false(check_rr_args(5, 5, 0.05, tst, 100, 0, TRUE)$ref.pvalue, info = tst)
  }
  expect_error(check_rr_args(5, 5, 0.05, 'Z-pool', 100, 0, NA), 'ref.pvalue')
  expect_error(check_rr_args(5, 5, 0.05, 'Z-pool', 100, 0, 1), 'ref.pvalue')
  expect_error(check_rr_args(5, 5, 0.05, 'Z-pool', 100, 0, c(TRUE, FALSE)), 'ref.pvalue')
})

test_that("every exported function with a grid takes ref.pvalue last, FALSE by default", {
  # devtools::load_all() exports the internal functions as well, so the exported functions
  # are identified by their common prefix
  ns <- asNamespace('bbssr')
  checked <- 0
  for (f in grep('^Binary', getNamespaceExports('bbssr'), value = TRUE)) {
    fun <- get(f, envir = ns)
    if (!is.function(fun) || !('n.grid' %in% names(formals(fun)))) next
    checked <- checked + 1
    a <- formals(fun)
    expect_identical(names(a)[length(a)], 'ref.pvalue', info = f)
    expect_identical(a$ref.pvalue, FALSE, info = f)
  }
  expect_gte(checked, 7)
})
