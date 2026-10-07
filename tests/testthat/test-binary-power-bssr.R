test_that("BinaryPowerBSSR returns a bbssr_powerbssr data frame", {
  res <- BinaryPowerBSSR(
    p = 0.45,
    Delta.A = 0.3, Delta.T = 0.3,
    N1 = 6, N2 = 6, omega = 0.5, r = 1,
    alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
  )
  expect_s3_class(res, 'bbssr_powerbssr')
  expect_named(res, c('p1', 'p2', 'p', 'power.BSSR', 'power.TRAD', 'E.N'))
  expect_equal(nrow(res), 1L)
  expect_true(all(res$power.BSSR >= 0 & res$power.BSSR <= 1))
  expect_true(all(res$power.TRAD >= 0 & res$power.TRAD <= 1))
})

test_that("the arguments of the weighted approach have been removed", {
  base <- list(p = 0.45, Delta.A = 0.3, Delta.T = 0.3, N1 = 6, N2 = 6,
               omega = 0.5, r = 1, alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
  for (dropped in list(list(weighted = TRUE), list(asmd.p1 = 0.6), list(asmd.p2 = 0.3))) {
    expect_error(do.call(BinaryPowerBSSR, c(base, dropped)), 'unused argument',
                 info = names(dropped))
  }
})

test_that("the remaining arguments are all used", {
  # A design is fully determined by Delta.A together with the initial sample sizes,
  # so no argument of the formals may be ignored by the body
  body.text <- paste(deparse(body(BinaryPowerBSSR)), collapse = ' ')
  for (arg in setdiff(names(formals(BinaryPowerBSSR)), '...')) {
    expect_match(body.text, paste0('(?<![\\w.])', arg, '(?![\\w.])'),
                 perl = TRUE, info = arg)
  }
})

test_that("BinaryPowerBSSR is vectorised over the pooled response probability", {
  res <- BinaryPowerBSSR(
    p = seq(0.3, 0.4, by = 0.05),
    Delta.A = 0.3, Delta.T = 0.3,
    N1 = 8, N2 = 8, omega = 0.5, r = 1,
    alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
  )
  expect_equal(nrow(res), 3L)
  expect_equal(res$p1 - res$p2, rep(0.3, 3))
})

test_that("the expected sample size is at least the interim sample size", {
  res <- BinaryPowerBSSR(
    p = 0.35,
    Delta.A = 0.3, Delta.T = 0.3,
    N1 = 10, N2 = 10, omega = 0.5, r = 1,
    alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
  )
  expect_gte(res$E.N, 10)
})

test_that("the restricted rule never enrols fewer patients than the unrestricted rule", {
  args <- list(p = c(0.3, 0.35),
               Delta.A = 0.3, Delta.T = 0.3, N1 = 12, N2 = 12, omega = 0.5, r = 1,
               alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
  unres <- do.call(BinaryPowerBSSR, c(args, list(restricted = FALSE)))
  res <- do.call(BinaryPowerBSSR, c(args, list(restricted = TRUE)))
  expect_true(all(res$E.N >= unres$E.N - 1e-9))
})

test_that("a true treatment effect of zero gives the type I error rate", {
  res <- BinaryPowerBSSR(
    p = c(0.3, 0.5),
    Delta.A = 0.3, Delta.T = 0,
    N1 = 10, N2 = 10, omega = 0.5, r = 1,
    alpha = 0.025, tar.power = 0.8, Test = 'Fisher'
  )
  expect_equal(res$p1, res$p2)
  expect_true(all(res$power.BSSR <= 0.05))
  expect_true(all(res$power.TRAD <= 0.025 + 1e-10))
})

test_that("scenarios outside the unit interval are dropped", {
  res <- BinaryPowerBSSR(
    p = c(0.2, 0.5, 0.95),
    Delta.A = 0.4, Delta.T = 0.4,
    N1 = 6, N2 = 6, omega = 0.5, r = 1,
    alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
  )
  expect_lt(nrow(res), 3L)
  expect_true(all(res$p1 <= 1 & res$p2 >= 0))
})

test_that("BinaryPowerBSSR accepts an interim fraction of one", {
  res <- BinaryPowerBSSR(
    p = 0.45,
    Delta.A = 0.3, Delta.T = 0.3,
    N1 = 6, N2 = 6, omega = 1, r = 1,
    alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
  )
  expect_true(is.finite(res$power.BSSR))
  expect_error(
    BinaryPowerBSSR(p = 0.45,
                    Delta.A = 0.3, Delta.T = 0.3, N1 = 6, N2 = 6, omega = 1.5, r = 1,
                    alpha = 0.025, tar.power = 0.8, Test = 'Chisq'),
    'omega'
  )
})

test_that("the allocation ratio is preserved by the interim and the final sample sizes", {
  for (r in c(1, 2, 3)) {
    res <- BinaryPowerBSSR(
      p = 0.35, Delta.A = 0.3, Delta.T = 0.3,
      N1 = ceiling(r * 11), N2 = 11, omega = 0.8, r = r,
      alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
    )
    # A whole allocation ratio must be reproduced exactly at the interim analysis
    expect_equal(attr(res, 'n1.interim'), as.integer(r * attr(res, 'n2.interim')),
                 info = sprintf('r = %g', r))
  }
})

test_that("a fractional allocation ratio keeps the interim sizes within one patient", {
  for (r in c(0.5, 1.5, 2.5)) {
    res <- BinaryPowerBSSR(
      p = 0.35, Delta.A = 0.3, Delta.T = 0.3,
      N1 = ceiling(r * 11), N2 = 11, omega = 0.8, r = r,
      alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
    )
    n1 <- attr(res, 'n1.interim')
    n2 <- attr(res, 'n2.interim')
    expect_equal(n1, as.integer(ceiling(r * n2)), info = sprintf('r = %g', r))
    expect_lt(abs(n1 - r * n2), 1)
  }
})

test_that("an initial size inconsistent with the allocation ratio is flagged", {
  expect_warning(
    BinaryPowerBSSR(p = 0.45, Delta.A = 0.3, Delta.T = 0.3,
                    N1 = 7, N2 = 6, omega = 0.5, r = 1,
                    alpha = 0.025, tar.power = 0.8, Test = 'Chisq'),
    'ratio r to 1'
  )
})

# Reference values for p = c(0.3, 0.45), Delta.A = 0.3, N2 = 10, omega = 0.5,
# alpha = 0.025 and tar.power = 0.8, computed with an independent Python implementation
# (numpy and scipy) that re-estimates the sample size by the same search and sums the
# joint binomial probabilities over the interim and the second-stage outcomes directly
bssr_ref <- list(
  list(Test = 'Chisq', alternative = 'greater', r = 1, Delta.T = 0.3,
       power = c(0.756413406488437, 0.750705544729915),
       E.N = c(65.3059377228023, 75.5787040704)),
  list(Test = 'Chisq', alternative = 'greater', r = 1, Delta.T = 0,
       power = c(0.0259608935849485, 0.0256688814937207),
       E.N = c(65.1643484471999, 74.8364847760476)),
  list(Test = 'Chisq', alternative = 'two.sided', r = 2, Delta.T = 0.3,
       power = c(0.77122764699632, 0.768379918422814),
       E.N = c(87.8982388035471, 104.98263754682)),
  list(Test = 'Chisq', alternative = 'two.sided', r = 2, Delta.T = 0,
       power = c(0.0249933395369282, 0.0263553983742971),
       E.N = c(88.2161329574932, 104.40801270423)),
  list(Test = 'Boschloo', alternative = 'greater', r = 2, Delta.T = 0.3,
       power = c(0.791254739800101, 0.77323324422443),
       E.N = c(74.9981051601408, 88.8779736124966)),
  list(Test = 'Boschloo', alternative = 'greater', r = 2, Delta.T = 0,
       power = c(0.021701985297022, 0.0236598768778865),
       E.N = c(75.3153533847524, 88.4155899013014)),
  list(Test = 'Z-pool', alternative = 'greater', r = 1, Delta.T = 0.3,
       power = c(0.761550296377532, 0.745360578459642),
       E.N = c(68.7037165753312, 78.3691243008)),
  list(Test = 'Z-pool', alternative = 'greater', r = 1, Delta.T = 0,
       power = c(0.0230522969203677, 0.0217694140545424),
       E.N = c(68.5263087743999, 77.6771381811188)),
  list(Test = 'Fisher', alternative = 'greater', r = 1, Delta.T = 0.3,
       power = c(0.777587216232051, 0.765350195503185),
       E.N = c(80.4580176190472, 89.5666384256)),
  list(Test = 'Fisher', alternative = 'greater', r = 1, Delta.T = 0,
       power = c(0.0142450507723455, 0.0158535596996759),
       E.N = c(80.4621801307999, 88.897260350634))
)

check_bssr_ref <- function(cases) {
  for (cs in cases) {
    res <- BinaryPowerBSSR(
      p = c(0.3, 0.45), Delta.A = 0.3, Delta.T = cs$Delta.T,
      N1 = ceiling(cs$r * 10), N2 = 10, omega = 0.5, r = cs$r,
      alpha = 0.025, tar.power = 0.8, Test = cs$Test, alternative = cs$alternative
    )
    info <- sprintf('%s, %s, r = %g, Delta.T = %g', cs$Test, cs$alternative, cs$r,
                    cs$Delta.T)
    expect_equal(res$power.BSSR, cs$power, tolerance = 1e-10, info = info)
    expect_equal(res$E.N, cs$E.N, tolerance = 1e-10, info = info)
  }
}

test_that("BinaryPowerBSSR agrees with an independent implementation (chi-squared)", {
  check_bssr_ref(Filter(function(cs) cs$Test == 'Chisq', bssr_ref))
})

test_that("BinaryPowerBSSR agrees with an independent implementation (exact tests)", {
  check_bssr_ref(Filter(function(cs) cs$Test != 'Chisq', bssr_ref))
})

test_that("the results do not depend on the reuse of rejection regions", {
  args <- list(p = c(0.3, 0.45), Delta.A = 0.3, Delta.T = 0.3, N1 = 10, N2 = 10,
               omega = 0.5, r = 1, alpha = 0.025, tar.power = 0.8, Test = 'Z-pool')
  clear_rr_cache()
  cold <- do.call(BinaryPowerBSSR, args)
  warm <- do.call(BinaryPowerBSSR, args)
  old <- options(bbssr.cache = FALSE)
  on.exit(options(old))
  off <- do.call(BinaryPowerBSSR, args)
  expect_identical(warm, cold)
  expect_identical(off, cold)
})

test_that("Friede and Kieser (2004), Table I is reproduced", {
  # The reference values were computed by the R package rbssr, an independent
  # implementation of the same design, and round to the published values
  fk <- BinaryPowerBSSR(
    p = 0.15, Delta.A = 0.2, Delta.T = 0.2, N1 = 51, N2 = 51, n.interim = c(20, 20),
    r = 1, alpha = 0.05, tar.power = 0.8, Test = 'Chisq', alternative = 'two.sided',
    ss.method = 'null.variance', rounding = 'friede-kieser', N.max = 204
  )
  expect_equal(fk$E.N, 99.0535867880608, tolerance = 1e-8)
  expect_equal(fk$power.BSSR, 0.778654816886751, tolerance = 1e-8)
  expect_equal(c(round(fk$E.N, 1), round(fk$power.BSSR, 3)), c(99.1, 0.779))
  fk3 <- BinaryPowerBSSR(
    p = 0.9, Delta.A = 0.2, Delta.T = 0.2, N1 = 71, N2 = 24, n.interim = c(30, 10),
    r = 3, alpha = 0.05, tar.power = 0.8, Test = 'Chisq', alternative = 'two.sided',
    ss.method = 'null.variance', rounding = 'friede-kieser', N.max = 190
  )
  expect_equal(fk3$E.N, 94.7067601478391, tolerance = 1e-8)
  expect_equal(fk3$power.BSSR, 0.693821084489667, tolerance = 1e-8)
  expect_equal(c(round(fk3$E.N, 1), round(fk3$power.BSSR, 3)), c(94.7, 0.694))
})

test_that("the lower-tail alternative mirrors the upper-tail one", {
  args <- list(p = c(0.3, 0.45), N1 = 10, N2 = 10, omega = 0.5, r = 1, alpha = 0.025,
               tar.power = 0.8, Test = 'Z-pool')
  up <- do.call(BinaryPowerBSSR, c(args, list(Delta.A = 0.3, Delta.T = 0.3)))
  down <- do.call(BinaryPowerBSSR, c(args, list(Delta.A = -0.3, Delta.T = -0.3,
                                                alternative = 'less')))
  expect_equal(down$power.BSSR, up$power.BSSR, tolerance = 1e-12)
  expect_equal(down$power.TRAD, up$power.TRAD, tolerance = 1e-12)
  expect_equal(down$E.N, up$E.N, tolerance = 1e-12)
})

test_that("n.interim gives the same design as the equivalent omega", {
  args <- list(p = c(0.3, 0.45), Delta.A = 0.3, Delta.T = 0.3, N1 = 10, N2 = 10, r = 1,
               alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
  a <- do.call(BinaryPowerBSSR, c(args, list(omega = 0.5)))
  b <- do.call(BinaryPowerBSSR, c(args, list(n.interim = c(5, 5))))
  expect_equal(b$power.BSSR, a$power.BSSR)
  expect_equal(b$E.N, a$E.N)
})

test_that("N.max bounds every final sample size", {
  res <- BinaryPowerBSSR(p = c(0.3, 0.5), Delta.A = 0.2, Delta.T = 0.2, N1 = 20, N2 = 20,
                         omega = 0.5, r = 1, alpha = 0.025, tar.power = 0.8,
                         Test = 'Chisq', ss.method = 'standard', N.max = 60)
  map <- attr(res, 'reestimation')
  expect_true(all(map$N1 + map$N2 <= 60))
  expect_true(any(map$N1 + map$N2 == 60))
  expect_true(all(res$E.N <= 60 + 1e-9))
})

test_that("a ratio effect splits the pooled proportion by the ratio", {
  res <- BinaryPowerBSSR(p = 0.3, Delta.A = 2, Delta.T = 2, N1 = 20, N2 = 20,
                         omega = 0.5, r = 1, alpha = 0.025, tar.power = 0.8,
                         Test = 'Chisq', effect = 'RR', ss.method = 'standard')
  expect_equal(res$p1 / res$p2, 2)
  expect_equal((res$p1 + res$p2) / 2, 0.3)
  map <- attr(res, 'reestimation')
  # No interim responder leaves both recovered rates at zero, so the plan is kept
  expect_equal(c(map$N1[1], map$N2[1]), c(20, 20))
})

test_that("scenarios outside the unit interval only by rounding error are evaluated", {
  # p[10] is 0.09999999999999999, so p2 = p - 0.1 is a rounding error below zero
  p <- seq(0.01, 0.5, by = 0.01)
  expect_no_warning(
    res <- BinaryPowerBSSR(p = p, Delta.A = 0.2, Delta.T = 0.2, N1 = 8, N2 = 8,
                           omega = 0.5, r = 1, alpha = 0.025, tar.power = 0.8,
                           Test = 'Chisq')
  )
  expect_equal(nrow(res), sum(p >= 0.1 - 1e-12))
  expect_true(all(res$p2 >= 0 & res$p1 <= 1))
  expect_true(all(is.finite(res$power.BSSR)))
})

test_that("BinaryPowerBSSR passes ref.pvalue to the analysis and the re-estimation", {
  run <- function(...) {
    BinaryPowerBSSR(p = c(0.3, 0.4), Delta.A = 0.3, Delta.T = 0.3, N1 = 14, N2 = 14,
                    omega = 0.5, r = 1, alpha = 0.025, tar.power = 0.8, Test = 'Z-pool',
                    ...)
  }
  seen <- ref_pvalue_calls(run(ref.pvalue = TRUE))
  expect_gt(nrow(seen), 1)
  expect_true(all(seen$ref.pvalue))
  expect_true(attr(run(ref.pvalue = TRUE), 'ref.pvalue'))
  # A conditional test used for the re-estimation has nothing to refine
  seen <- ref_pvalue_calls(run(ss.Test = 'Chisq', ref.pvalue = TRUE))
  expect_true(all(seen$ref.pvalue == (seen$Test == 'Z-pool')))
  expect_true(any(seen$Test == 'Chisq'))
  seen <- ref_pvalue_calls(run())
  expect_gt(nrow(seen), 0)
  expect_false(any(seen$ref.pvalue))
})

test_that("BinaryPowerBSSR handles a non-inferiority margin", {
  # Interim analysis of 30 per group, margin 0.2, re-estimation by the formula of
  # Farrington and Manning with each group rounded to the nearest whole number. Reference
  # values from tools/reference/reference_values.py
  run <- function(p, Delta.T) {
    BinaryPowerBSSR(p = p, Delta.A = 0, Delta.T = Delta.T, N1 = 54, N2 = 54,
                    n.interim = c(30, 30), r = 1, alpha = 0.025, tar.power = 0.8,
                    Test = 'Farrington-Manning', ss.method = 'standard',
                    rounding = 'nearest', margin = 0.2)
  }
  seen <- ref_pvalue_calls(res <- run(0.4, 0))
  expect_gt(nrow(seen), 0)
  expect_true(all(seen$margin == 0.2))
  expect_equal(attr(res, 'margin'), 0.2)
  expect_equal(c(res$power.BSSR, res$E.N), c(0.792565321803051, 181.159310064003),
               tolerance = 1e-10)
  # On the null boundary the power is the type I error rate
  res <- run(0.5, -0.2)
  expect_equal(c(res$p1, res$p2), c(0.4, 0.6))
  expect_equal(c(res$power.BSSR, res$E.N), c(0.0245342677114993, 187.723544195855),
               tolerance = 1e-10)
  expect_error(BinaryPowerBSSR(p = 0.4, Delta.A = -0.2, Delta.T = 0, N1 = 54, N2 = 54,
                               n.interim = c(30, 30), r = 1, alpha = 0.025,
                               tar.power = 0.8, Test = 'Farrington-Manning',
                               ss.method = 'standard', margin = 0.2), 'exceed -margin')
})

test_that("BinaryPowerBSSR passes tsmethod to the analysis and the re-estimation", {
  seen <- ref_pvalue_calls(BinaryPowerBSSR(p = 0.4, Delta.A = 0.3, Delta.T = 0.3, N1 = 16,
                                           N2 = 8, omega = 0.5, r = 2, alpha = 0.05,
                                           tar.power = 0.8, Test = 'Fisher',
                                           alternative = 'two.sided',
                                           tsmethod = 'blaker'))
  expect_gt(nrow(seen), 1)
  expect_true(all(seen$tsmethod == 'blaker'))
})

test_that("the exact re-estimation uses the search of BinarySampleSize", {
  args <- list(p = 0.4, Delta.A = 0.3, Delta.T = 0.3, N1 = 20, N2 = 20,
               n.interim = c(10, 10), r = 1, alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
  map <- lapply(c('crossing', 'smallest', 'stable'), function(s) {
    attr(do.call(BinaryPowerBSSR, c(args, list(search = s))), 'reestimation')
  })
  # Reference values from tools/reference/reference_values.py: the re-estimated size of
  # group 2 for every pooled interim count s = 0, ..., 20. The searches 'crossing' and
  # 'smallest' differ at s = 10, and 'stable' agrees with 'crossing' here
  expect_equal(map[[1]]$N2.re, c(36, 27, 21, 18, 24, 31, 34, 39, 40, 41, 44, 41, 40, 39, 34,
                                 31, 24, 18, 21, 27, 36))
  expect_equal(map[[2]]$N2.re, c(36, 27, 21, 18, 24, 31, 34, 39, 40, 41, 41, 41, 40, 39, 34,
                                 31, 24, 18, 21, 27, 36))
  expect_equal(map[[3]]$N2.re, map[[1]]$N2.re)
  expect_equal(map[[3]]$N2.limit, c(98, 85, 77, 72, 77, 82, 86, 89, 91, 92, 93, 92, 91, 89,
                                    86, 82, 77, 72, 77, 85, 98))
  expect_true(all(is.na(map[[1]]$N2.limit) & is.na(map[[2]]$N2.limit)))
})

test_that("BinaryPowerBSSR passes search and search.limit to the exact re-estimation", {
  seen <- search_calls(res <- BinaryPowerBSSR(
    p = 0.4, Delta.A = 0.3, Delta.T = 0.3, N1 = 8, N2 = 8, n.interim = c(4, 4), r = 1,
    alpha = 0.025, tar.power = 0.8, Test = 'Chisq', search = 'stable',
    search.limit = c(3, 10)
  ))
  expect_gt(nrow(seen), 0)
  expect_true(all(seen$search == 'stable' & seen$a == 3 & seen$b == 10))
  expect_equal(attr(res, 'search'), 'stable')
  expect_false(anyNA(attr(res, 'reestimation')$N2.limit))
})
