ptie <- function(...) {
  args <- list(Delta.A = 0.25, N = 30, n.interim = 10, alpha = 0.05, tar.power = 0.8,
               Test = 'McNemar', ss.method = 'standard')
  new <- list(...)
  args[names(new)] <- new
  do.call(PairedTypeIErrorBSSR, args)
}

test_that("PairedTypeIErrorBSSR agrees with the reference values on a grid", {
  tie <- ptie(psi = c(0.1, 0.4, 0.7), maximize = 'grid')
  expect_equal(c(tie$TIE.BSSR, tie$TIE.TRAD),
               c(0.0360003775113387, 0.0500079876956218, 0.0494240699772004,
                 0.0460091085323312, 0.0516879128204756, 0.0485609956965745),
               tolerance = 1e-10)
  expect_equal(tie$p10, tie$p01)
  tie <- ptie(Delta.A = 0, Test = 'Paired-Score', margin = 0.2, psi = c(0.28, 0.52, 0.76),
              maximize = 'grid')
  expect_equal(c(tie$TIE.BSSR, tie$TIE.TRAD),
               c(0.0763047647322819, 0.0535085909814677, 0.0519467569184633,
                 0.0528469342611255, 0.0503274423521533, 0.0533225288037961),
               tolerance = 1e-10)
  expect_equal(tie$p10 - tie$p01, rep(-0.2, 3))
})

test_that("PairedTypeIErrorBSSR certifies the largest rates", {
  m <- attr(ptie(), 'max')
  expect_equal(m$TIE, c(0.0525269607489134, 0.0558117064555376), tolerance = 1e-9)
  expect_true(all(m$bound >= m$TIE & m$bound - m$TIE <= 1e-12 + 1e-15))
  m <- attr(ptie(Delta.A = 0, Test = 'Paired-Score', margin = 0.2), 'max')
  expect_equal(m$TIE, c(0.100021274533563, 0.0550269955431753), tolerance = 1e-9)
  expect_true(all(m$psi >= 0.2))
})

test_that("PairedTypeIErrorBSSR drops the values of psi off the boundary", {
  tie <- ptie(Delta.A = 0, Test = 'Paired-Score', margin = 0.2,
              psi = c(0, 0.1, 0.3, 0.6), maximize = 'grid')
  expect_equal(tie$psi, c(0.3, 0.6))
  expect_error(ptie(Delta.A = 0, Test = 'Paired-Score', margin = 0.2, psi = 0.1),
               'no value of psi')
  expect_error(ptie(psi = 1.5), 'psi must be')
})

test_that("PairedTypeIErrorBSSR evaluates both boundaries of a test of equivalence", {
  tie <- ptie(Delta.A = 0, Test = 'Paired-Score', alternative = 'equivalence',
              margin = 0.3, psi = c(0.4, 0.6, 0.8), maximize = 'grid')
  expect_equal(tie$boundary, rep(c('lower', 'upper'), each = 3))
  expect_equal(tie$p10 - tie$p01, rep(c(-0.3, 0.3), each = 3))
  # The two boundaries are mirror images, so the lower one is reported
  m <- attr(tie, 'max')
  expect_identical(m$boundary, c('lower', 'lower'))
})

test_that("the exact conditional McNemar test keeps its level after re-estimation", {
  m <- attr(ptie(Test = 'McNemar-exact', alternative = 'two.sided',
                 ss.method = 'conditional.power', N.max = 45), 'max')
  expect_true(all(m$bound <= 0.05))
})
