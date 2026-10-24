# Reference values are computed with an independent Python implementation
# (tools/reference/reference_strat.py, cond_reject) that sums the hypergeometric
# distributions of the group 1 responders directly. The design has the planned sizes 6 per
# group and stratum, the interim sizes 3, Delta.A = 0.35 on the risk difference, the level
# 0.05, re-estimation by the normal approximation and N.max = 36
crp <- function(...) {
  args <- list(Delta.A = 0.35, N1 = c(6, 6), N2 = c(6, 6), omega = 0.5, r = 1,
               alpha = 0.05, tar.power = 0.8, Test = 'CMH', ss.method = 'standard',
               N.max = 36, theta = rbind(c(0.3, 0.5), c(0.5, 0.5)))
  new <- list(...)
  args[names(new)] <- new
  do.call(StratCondRejectBSSR, args)
}

test_that("StratCondRejectBSSR agrees with an independent implementation", {
  ref <- list(
    CMH = list(sums = c(5.04246914881655, 0.0321590909090909, 4.85014382410244),
               top = c(0.13, 0.121212121212121),
               loc = c(2, 0, 3, 0),
               exceed = c(1497, 2092),
               tie = c(0.0555498675564424, 0.0564531616983004),
               cell = c(0.0631218753516136, 0.0595825994014741)),
    'Strat-Fisher' = list(sums = c(1.85932951968666, 0.0115909090909091, 2.17966362655884),
                          top = c(0.0553617576507187, 0.0497077580633126),
                          loc = c(3, 3, 6, 6),
                          exceed = c(61, 0),
                          tie = c(0.026349572048717, 0.0316584061947651),
                          cell = c(0.0320588949234571, 0.0357064372301102)),
    'Unstrat-Chisq' = list(sums = c(4.33050910894661, 0.0225, 4.73128420835368),
                           top = c(0.125, 0.117647058823529),
                           loc = c(1, 0, 1, 1),
                           exceed = c(1849, 2864),
                           tie = c(0.0583615622080922, 0.0658128868963104),
                           cell = c(0.0653565634470866, 0.0706095910072328)))
  for (Test in names(ref)) {
    res <- crp(Test = Test)
    r <- ref[[Test]]
    expect_equal(c(sum(res$CRP.max), sum(res$CRP.min), sum(res$CRP.total.max)), r$sums,
                 tolerance = 1e-10)
    m <- attr(res, 'max')
    expect_equal(m$value, r$top, tolerance = 1e-10)
    expect_equal(c(m$s1[1], m$s2[1], m$w1[1], m$w2[1]), r$loc)
    expect_equal(unname(attr(res, 'exceed')), r$exceed)
    expect_equal(attr(res, 'TIE')$TIE, r$tie, tolerance = 1e-10)
    b <- attr(res, 'by.s')
    expect_equal(b$TIE.s[b$s1 == 2 & b$s2 == 3], r$cell, tolerance = 1e-10)
  }
})

test_that("StratCondRejectBSSR decomposes the type I error rate of StratTypeIErrorBSSR", {
  th <- rbind(c(0.3, 0.5), c(0.5, 0.5), c(0.85, 0.2))
  res <- crp(theta = th)
  tie <- StratTypeIErrorBSSR(Delta.A = 0.35, N1 = c(6, 6), N2 = c(6, 6), omega = 0.5, r = 1,
                             alpha = 0.05, tar.power = 0.8, Test = 'CMH',
                             ss.method = 'standard', N.max = 36, theta = th,
                             maximize = 'grid')
  expect_equal(attr(res, 'TIE')$TIE, tie$TIE.BSSR, tolerance = 1e-12)
  b <- attr(res, 'by.s')
  expect_identical(nrow(b), 3L * nrow(res))
  # The probabilities of the interim outcomes sum to 1 at every point
  expect_equal(as.vector(tapply(b$prob.s, b$theta.S1, sum)), c(1, 1, 1))
  expect_equal(b$contribution, b$prob.s * b$TIE.s)
})

test_that("StratCondRejectBSSR lists the interim outcomes and the full listing", {
  res <- crp(Test = 'Strat-Fisher', full = TRUE, theta = NULL)
  expect_s3_class(res, 'bbssr_scrp')
  expect_named(res, c('s1', 's2', 'N1.S1', 'N1.S2', 'N2.S1', 'N2.S2', 'CRP.max', 'CRP.min',
                      'exceed', 'CRP.total.max'))
  expect_identical(nrow(res), 7L * 7L)
  expect_identical(res$s1[1:8], c(0:6, 0L))
  expect_null(attr(res, 'by.s'))
  map <- attr(res, 'reestimation')
  expect_identical(res$N1.S2, map$N1.S2)
  f <- attr(res, 'full')
  expect_identical(nrow(f), as.integer(sum((res$N1.S1 + res$N2.S1 - 5) *
                                             (res$N1.S2 + res$N2.S2 - 5))))
  expect_identical(order(f$s2, f$s1, f$w2, f$w1), seq_len(nrow(f)))
  key <- paste(f$s1, f$s2)
  expect_equal(as.vector(tapply(f$CRP, key, max)[paste(res$s1, res$s2)]), res$CRP.max)
  expect_equal(as.vector(tapply(f$CRP, key, min)[paste(res$s1, res$s2)]), res$CRP.min)
  expect_equal(sum(f$CRP - 0.05 > 1.4901161193847656e-08), attr(res, 'exceed')[['CRP']])
  # The stratified exact test keeps the conditional size given the totals at the level
  expect_lte(max(f$CRP.total), 0.05)
  expect_identical(attr(res, 'exceed')[['CRP.total']], 0L)
})

test_that("StratCondRejectBSSR checks its arguments", {
  expect_error(crp(N1 = c(6, 6, 6), N2 = c(6, 6, 6)), 'two strata')
  expect_error(crp(full = NA), 'full must be')
  expect_error(crp(theta = c(0.2, 1.2)))
  expect_error(crp(N1 = c(60, 60), N2 = c(60, 60), N.max = 480, full = TRUE, theta = NULL,
                   ss.method = 'standard'), 'more than 1e6')
})

