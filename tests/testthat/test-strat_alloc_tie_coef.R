test_that("strat_alloc_tie_coef gives the exact type I error rate of a random allocation", {
  # The polynomial at a point equals the power of StratPower with the same allocation
  # under the null hypothesis, computed by another route
  pts <- rbind(c(0.3, 0.6), c(0.75, 0.2), c(0.5, 0.5))
  for (alloc in c('stratum', 'group', 'both')) {
    for (Test in c('CMH', 'Strat-Fisher', 'Unstrat-Chisq')) {
      al <- check_strat_alloc(alloc, 16, c(0.5, 0.5), 1, 4)
      res <- strat_alloc_tie_coef(16, al, alloc, 0.05, Test, 'greater', 'minlike')
      pw <- StratPower(pts, pts, alpha = 0.05, Test = Test, alloc = alloc, N = 16,
                       frac = c(0.5, 0.5))$Power
      expect_equal(strat_bern_eval(res$coef, pts[, 1], pts[, 2]), pw, tolerance = 1e-10,
                   info = paste(alloc, Test))
    }
  }
  # The degrees are the largest stratum sizes, the fixed sizes with alloc = 'group'
  al <- check_strat_alloc('group', 17, c(0.4, 0.6), 1, 4)
  res <- strat_alloc_tie_coef(17, al, 'group', 0.05, 'CMH', 'greater', 'minlike')
  expect_identical(dim(res$coef), c(7L, 12L))
})
