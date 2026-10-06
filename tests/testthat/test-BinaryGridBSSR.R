# Reference values from tools/reference/reference_values.py (block_grid), an independent
# Python implementation that re-estimates the sample size by the same rules and sums the
# binomial probabilities over the interim and the second-stage outcomes directly
test_that("BinaryGridBSSR agrees with an independent implementation", {
  design <- data.frame(Test = c('Chisq', 'Fisher'))
  res <- BinaryGridBSSR(design, p = c(0.3, 0.45), Delta.A = 0.3, N1 = 10, N2 = 10,
                        omega = 0.5, r = 1, alpha = 0.025, tar.power = 0.8)
  expect_s3_class(res, 'bbssr_grid')
  expect_equal(res$design, c(1, 1, 2, 2))
  expect_equal(res$Test, c('Chisq', 'Chisq', 'Fisher', 'Fisher'))
  expect_equal(res$Delta.T, rep(0.3, 4))
  expect_equal(res$power.BSSR,
               c(0.756413406488437, 0.750705544729915, 0.777587216232051,
                 0.765350195503184), tolerance = 1e-10)
  expect_equal(res$E.N,
               c(65.3059377228023, 75.5787040704, 80.4580176190472, 89.5666384256),
               tolerance = 1e-10)
  expect_equal(res$SD.N,
               c(15.4108334577559, 12.1779563802348, 14.8447550877925,
                 11.3993961793173), tolerance = 1e-10)
})

test_that("each design gives the result of BinaryPowerBSSR and its summary", {
  # The second design gives the interim sizes instead of omega
  design <- data.frame(Test = c('Chisq', 'Z-pool'), omega = c(0.5, NA),
                       n1.interim = c(NA, 6), n2.interim = c(NA, 6))
  p <- c(0.25, 0.4, 0.55)
  res <- BinaryGridBSSR(design, p = p, Delta.A = 0.3, N1 = 14, N2 = 14, r = 1,
                        alpha = 0.025, tar.power = 0.8, ss.method = 'standard')
  direct <- list(
    BinaryPowerBSSR(p = p, Delta.A = 0.3, Delta.T = 0.3, N1 = 14, N2 = 14, omega = 0.5,
                    r = 1, alpha = 0.025, tar.power = 0.8, Test = 'Chisq',
                    ss.method = 'standard'),
    BinaryPowerBSSR(p = p, Delta.A = 0.3, Delta.T = 0.3, N1 = 14, N2 = 14,
                    n.interim = c(6, 6), r = 1, alpha = 0.025, tar.power = 0.8,
                    Test = 'Z-pool', ss.method = 'standard')
  )
  for (i in 1:2) {
    k <- res$design == i
    d <- direct[[i]]
    s <- summary(d)
    expect_equal(sum(k), 3)
    expect_identical(res$p[k], d$p)
    expect_identical(res$p1[k], d$p1)
    expect_identical(res$power.BSSR[k], d$power.BSSR)
    expect_identical(res$power.TRAD[k], d$power.TRAD)
    expect_identical(res$E.N[k], d$E.N)
    expect_identical(res$SD.N[k], s$SD.N)
    expect_identical(res$N.q50[k], s$N.q50)
    expect_identical(res$P.N.max[k], s$P.N.max)
    expect_identical(unique(res$n1.interim[k]), attr(d, 'n1.interim'))
  }
  expect_equal(nrow(attr(res, 'designs')), 2)
  # A factor column, as made by expand.grid(), gives the same results
  fac <- BinaryGridBSSR(expand.grid(Test = c('Chisq', 'Fisher')), p = 0.4, Delta.A = 0.3,
                        N1 = 14, N2 = 14, omega = 0.5, r = 1, alpha = 0.025,
                        tar.power = 0.8, ss.method = 'standard')
  chr <- BinaryGridBSSR(data.frame(Test = c('Chisq', 'Fisher')), p = 0.4, Delta.A = 0.3,
                        N1 = 14, N2 = 14, omega = 0.5, r = 1, alpha = 0.025,
                        tar.power = 0.8, ss.method = 'standard')
  expect_identical(fac$power.BSSR, chr$power.BSSR)
})

test_that("the initial sample sizes follow from p.plan", {
  design <- data.frame(ss.method = c('exact', 'standard'), r = c(1, 2), p.plan = 0.4)
  res <- BinaryGridBSSR(design, p = 0.4, Delta.A = 0.3, omega = 0.5, alpha = 0.025,
                        tar.power = 0.8, Test = 'Chisq')
  d <- attr(res, 'designs')
  # 40 per group from the exact search, and 60 and 30 from the normal approximation
  expect_equal(c(d$N1[1], d$N2[1]), c(40, 40))
  expect_equal(c(d$N1[2], d$N2[2]), c(60, 30))
  ss <- BinarySampleSize(0.55, 0.25, 1, 0.025, 0.8, 'Chisq')
  expect_identical(c(d$N1[1], d$N2[1]), c(ss$N1, ss$N2))
  ss <- BinarySampleSize(0.5, 0.2, 2, 0.025, 0.8, 'Chisq', method = 'standard')
  expect_identical(c(d$N1[2], d$N2[2]), c(ss$N1, ss$N2))
  expect_equal(res$power.TRAD[1], 0.801503911679796, tolerance = 1e-10)
})

test_that("Delta.T is taken from the argument, a column or the assumed effect", {
  run <- function(design, ...) {
    BinaryGridBSSR(design, p = 0.4, Delta.A = 0.3, N1 = 12, N2 = 12, omega = 0.5, r = 1,
                   tar.power = 0.8, Test = 'Chisq', ss.method = 'standard', ...)
  }
  a <- run(data.frame(alpha = 0.025))
  b <- run(data.frame(alpha = 0.025), Delta.T = 0)
  d <- run(data.frame(alpha = 0.025, Delta.T = 0))
  expect_equal(a$Delta.T, 0.3)
  expect_equal(b$Delta.T, 0)
  direct <- BinaryPowerBSSR(p = 0.4, Delta.A = 0.3, Delta.T = 0, N1 = 12, N2 = 12,
                            omega = 0.5, r = 1, alpha = 0.025, tar.power = 0.8,
                            Test = 'Chisq', ss.method = 'standard')
  expect_identical(b$power.BSSR, direct$power.BSSR)
  expect_identical(d$power.BSSR, direct$power.BSSR)
  expect_error(run(data.frame(alpha = 0.025, Delta.T = 0), Delta.T = 0),
               'either as the argument')
})

test_that("type1 adds the largest type I error rates of BinaryTypeIErrorBSSR", {
  design <- data.frame(Delta.A = c(0.3, 0.4))
  res <- BinaryGridBSSR(design, p = 0.4, N1 = 12, N2 = 12, omega = 0.5, r = 1,
                        alpha = 0.025, tar.power = 0.8, Test = 'Chisq',
                        ss.method = 'standard', type1 = TRUE)
  d <- attr(res, 'designs')
  expect_named(d, c('design', 'Delta.A', 'N1', 'N2', 'n1.interim', 'n2.interim',
                    'Delta.T', 'TIE.BSSR', 'theta.BSSR', 'TIE.TRAD', 'theta.TRAD'))
  for (i in 1:2) {
    tie <- BinaryTypeIErrorBSSR(Delta.A = design$Delta.A[i], N1 = 12, N2 = 12,
                                omega = 0.5, r = 1, alpha = 0.025, tar.power = 0.8,
                                Test = 'Chisq', ss.method = 'standard')
    m <- attr(tie, 'max')
    expect_identical(c(d$TIE.BSSR[i], d$TIE.TRAD[i]), m$TIE)
    expect_identical(c(d$theta.BSSR[i], d$theta.TRAD[i]), m$theta)
  }
  expect_true(attr(res, 'type1'))
})

test_that("BinaryGridBSSR reports the design whose arguments are wrong", {
  run <- function(design, ...) {
    BinaryGridBSSR(design, p = 0.4, Delta.A = 0.3, r = 1, alpha = 0.025, tar.power = 0.8,
                   Test = 'Chisq', ss.method = 'standard', ...)
  }
  expect_error(run(data.frame(omega = 0.5, N1 = 12, N2 = 12, foo = 1)),
               'unknown argument\\(s\\): foo')
  expect_error(run(data.frame(omega = 0.5, N1 = 12, N2 = 12, r = 1)), 'given both')
  expect_error(run(data.frame(omega = 0.5, N1 = c(12, NA), N2 = c(12, NA))),
               'design 2: give either N1 and N2 or p.plan')
  expect_error(run(data.frame(omega = 0.5, N1 = 12, N2 = 12, p.plan = 0.4)),
               'design 1: give either N1 and N2 or p.plan, not both')
  expect_error(run(data.frame(N1 = 12, N2 = 12, n1.interim = 6)),
               'design 1: give both n1.interim and n2.interim')
  expect_error(run(data.frame(omega = 0.5, N1 = 12)), 'give both N1 and N2')
  expect_error(BinaryGridBSSR(data.frame(Test = 'Chisq'), p = 0.4, Delta.A = 0.3, N1 = 12,
                              N2 = 12, omega = 0.5, alpha = 0.025, tar.power = 0.8),
               'design 1: r must be given')
  expect_error(BinaryGridBSSR(data.frame(), p = 0.4), 'at least one row')
  expect_error(run(data.frame(omega = 0.5, N1 = 12, N2 = 12), Delta.T = c(0, 0.3)),
               'design 1: Delta.T must be a single number')
  # Every value of p puts p2 below 0 for the second design
  expect_error(BinaryGridBSSR(data.frame(Delta.A = c(0.3, 0.5)), p = c(0.15, 0.2),
                              N1 = 12, N2 = 12, omega = 0.5, r = 1, alpha = 0.025,
                              tar.power = 0.8, Test = 'Chisq', ss.method = 'standard'),
               'design 2: no scenario has both p1 and p2 inside the unit interval')
})

test_that("a warning names the design that issued it", {
  expect_warning(
    BinaryGridBSSR(data.frame(N1 = 13), p = 0.4, Delta.A = 0.3, N2 = 12, omega = 0.5,
                   r = 1, alpha = 0.025, tar.power = 0.8, Test = 'Chisq',
                   ss.method = 'standard'),
    'design 1: N1 differs from ceiling'
  )
})

test_that("verbose reports each design", {
  expect_message(
    BinaryGridBSSR(data.frame(Delta.A = 0.3), p = 0.4, N1 = 12, N2 = 12, omega = 0.5,
                   r = 1, alpha = 0.025, tar.power = 0.8, Test = 'Chisq',
                   ss.method = 'standard', verbose = TRUE),
    'design 1 of 1'
  )
})

test_that("BinaryGridBSSR passes tsmethod given as a column of the design", {
  design <- data.frame(tsmethod = c('minlike', 'blaker'), stringsAsFactors = FALSE)
  seen <- ref_pvalue_calls(g <- BinaryGridBSSR(design, p = 0.4, Delta.A = 0.3, N1 = 16,
                                               N2 = 8, omega = 0.5, r = 2, alpha = 0.05,
                                               tar.power = 0.8, Test = 'Fisher',
                                               alternative = 'two.sided',
                                               ss.method = 'standard'))
  expect_setequal(unique(seen$tsmethod), c('minlike', 'blaker'))
  expect_equal(g$tsmethod, c('minlike', 'blaker'))
})
