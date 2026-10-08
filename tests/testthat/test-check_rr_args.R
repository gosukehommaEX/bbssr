# check_rr_args() with every argument named and a default for each, so that each
# expectation varies one argument only
args_of <- function(N1 = 5, N2 = 5, alpha = 0.05, Test = 'Chisq', n.grid = 100,
                    bb.gamma = 0, ref.pvalue = FALSE, alternative = 'greater', margin = 0,
                    margin.scale = 'RD') {
  check_rr_args(N1, N2, alpha, Test, n.grid, bb.gamma, ref.pvalue, alternative, margin,
                margin.scale)
}

test_that("check_rr_args returns the arguments in computational form", {
  a <- check_rr_args(10, 8, 0.025, 'Bosch', 50, 0, TRUE, 'greater', 0, 'RD')
  expect_identical(a, list(N1 = 10L, N2 = 8L, n.grid = 50L, Test = 'Boschloo',
                           ref.pvalue = TRUE, margin = 0, margin.scale = 'RD'))
  expect_identical(args_of(Test = 'Farr', margin = 0.1)$Test, 'Farrington-Manning')
  expect_identical(args_of(Test = 'Black', margin = -0.2)$margin, -0.2)
})

test_that("check_rr_args rejects invalid arguments", {
  expect_error(args_of(N1 = c(5, 6)), 'single value')
  expect_error(args_of(N1 = 0), 'positive integers')
  expect_error(args_of(N2 = 4.5), 'positive integers')
  expect_error(args_of(alpha = 0), 'alpha')
  expect_error(args_of(n.grid = 1), 'n.grid')
  expect_error(args_of(Test = 'Z-pool', bb.gamma = 0.05), 'bb.gamma')
  expect_error(args_of(Test = 'Student'), 'should be one of')
})

test_that("check_rr_args warns when bb.gamma is given for a conditional test", {
  expect_warning(args_of(Test = 'Fisher', bb.gamma = 0.001), 'bb.gamma')
  expect_no_warning(args_of(Test = 'Z-pool', bb.gamma = 0.001))
})

test_that("check_rr_args validates ref.pvalue and keeps it for unconditional tests", {
  expect_true(args_of(Test = 'Z-pool', ref.pvalue = TRUE)$ref.pvalue)
  expect_true(args_of(Test = 'Boschloo', ref.pvalue = TRUE)$ref.pvalue)
  for (tst in c('Chisq', 'Fisher', 'Fisher-midP', 'Blackwelder', 'Farrington-Manning')) {
    expect_false(args_of(Test = tst, ref.pvalue = TRUE)$ref.pvalue, info = tst)
  }
  expect_error(args_of(Test = 'Z-pool', ref.pvalue = NA), 'ref.pvalue')
  expect_error(args_of(Test = 'Z-pool', ref.pvalue = 1), 'ref.pvalue')
  expect_error(args_of(Test = 'Z-pool', ref.pvalue = c(TRUE, FALSE)), 'ref.pvalue')
})

test_that("check_rr_args validates the margin", {
  expect_error(args_of(Test = 'Blackwelder', margin = NA), 'margin')
  expect_error(args_of(Test = 'Blackwelder', margin = 1), 'margin')
  expect_error(args_of(Test = 'Blackwelder', margin = c(0.1, 0.2)), 'margin')
  expect_error(args_of(Test = 'Blackwelder', margin = '0.1'), 'margin')
  for (tst in c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')) {
    expect_error(args_of(Test = tst, margin = 0.1), 'Blackwelder', info = tst)
  }
  expect_error(args_of(Test = 'Blackwelder', margin = 0.1, alternative = 'two.sided'),
               'one-sided')
  expect_no_error(args_of(Test = 'Blackwelder', alternative = 'two.sided'))
})

test_that("check_rr_args validates a margin on the scale of the risk ratio", {
  a <- args_of(Test = 'Farr', margin = 1.25, margin.scale = 'RR', alternative = 'less')
  expect_identical(a$margin, 1.25)
  expect_identical(a$margin.scale, 'RR')
  expect_identical(args_of(Test = 'Black', margin = 1, margin.scale = 'RR')$margin, 1)
  for (m in list(0, -0.2, NA, Inf, '0.8', c(0.8, 0.9))) {
    expect_error(args_of(Test = 'Blackwelder', margin = m, margin.scale = 'RR'),
                 'positive value', info = paste(m, collapse = ' '))
  }
  for (tst in c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')) {
    expect_error(args_of(Test = tst, margin = 0.8, margin.scale = 'RR'), 'Blackwelder',
                 info = tst)
  }
  expect_error(args_of(Test = 'Blackwelder', margin = 0.8, margin.scale = 'RR',
                       alternative = 'two.sided'), 'one-sided')
  expect_error(args_of(Test = 'Blackwelder', margin = 0.8, margin.scale = 'OR'),
               'margin.scale')
  expect_error(args_of(Test = 'Blackwelder', margin = 0.8, margin.scale = c('RD', 'RR')),
               'margin.scale')
})

test_that("every exported function with a grid ends with margin, ref.pvalue, margin.scale", {
  # devtools::load_all() exports the internal functions as well, so the exported functions
  # are identified by their common prefix
  ns <- asNamespace('bbssr')
  checked <- 0
  with.scale <- 0
  for (f in grep('^Binary', getNamespaceExports('bbssr'), value = TRUE)) {
    fun <- get(f, envir = ns)
    if (!is.function(fun) || !('n.grid' %in% names(formals(fun)))) next
    checked <- checked + 1
    a <- formals(fun)
    # BinaryCondRejectBSSR covers superiority only and has no margin.scale
    if ('margin.scale' %in% names(a)) {
      with.scale <- with.scale + 1
      expect_identical(names(a)[length(a)], 'margin.scale', info = f)
      expect_identical(eval(a$margin.scale), c('RD', 'RR'), info = f)
      a <- a[-length(a)]
    }
    expect_identical(names(a)[length(a)], 'ref.pvalue', info = f)
    expect_identical(a$ref.pvalue, FALSE, info = f)
    # The margin comes just before ref.pvalue and is 0 by default
    expect_identical(names(a)[length(a) - 1L], 'margin', info = f)
    expect_identical(a$margin, 0, info = f)
  }
  expect_gte(checked, 8)
  # BinaryRR, BinaryPower, BinarySampleSize, BinaryPowerBSSR, BinaryBSSR,
  # BinaryTypeIErrorBSSR and BinaryAlphaAdjBSSR
  expect_equal(with.scale, 7)
})

test_that("every exported function with tsmethod offers the three conventions", {
  ns <- asNamespace('bbssr')
  checked <- 0
  for (f in grep('^Binary', getNamespaceExports('bbssr'), value = TRUE)) {
    fun <- get(f, envir = ns)
    if (!is.function(fun) || !('tsmethod' %in% names(formals(fun)))) next
    checked <- checked + 1
    expect_identical(eval(formals(fun)$tsmethod), c('minlike', 'central', 'blaker'),
                     info = f)
  }
  expect_equal(checked, 8)
})

test_that("every exported function with search offers the three searches", {
  ns <- asNamespace('bbssr')
  checked <- 0
  for (f in grep('^Binary', getNamespaceExports('bbssr'), value = TRUE)) {
    fun <- get(f, envir = ns)
    if (!is.function(fun) || !('search' %in% names(formals(fun)))) next
    checked <- checked + 1
    expect_identical(eval(formals(fun)$search), c('crossing', 'smallest', 'stable'),
                     info = f)
    expect_identical(eval(formals(fun)$search.limit), c(2, 50), info = f)
  }
  # BinarySampleSize, BinaryPowerBSSR, BinaryBSSR, BinaryTypeIErrorBSSR,
  # BinaryAlphaAdjBSSR and BinaryCondRejectBSSR
  expect_equal(checked, 6)
})
