adj <- function(...) {
  args <- list(Delta.A = 0.35, N1 = c(6, 6), N2 = c(6, 6), omega = 0.5, r = 1,
               alpha = 0.1, tar.power = 0.8, Test = 'CMH', ss.method = 'standard',
               theta = seq(0.1, 0.9, by = 0.2), maximize = 'grid')
  new <- list(...)
  args[names(new)] <- new
  do.call(StratAlphaAdjBSSR, args)
}

test_that("StratAlphaAdjBSSR agrees with the reference levels on the grid", {
  res <- adj()
  expect_equal(res$max.TIE, c(0.108801399881486, 0.109288192949331), tolerance = 1e-10)
  expect_equal(res$alpha.adj, c(0.0947488, 0.0786496), tolerance = 1e-12)
  expect_equal(res$max.TIE.adj, c(0.0954459061613307, 0.0905266881854914),
               tolerance = 1e-10)
  expect_s3_class(res, 'bbssr_alphaadj')
  expect_named(res, c('Design', 'alpha', 'max.TIE', 'alpha.adj', 'max.TIE.adj',
                      'theta.S1.adj', 'theta.S2.adj', 'max.TIE.bound',
                      'max.TIE.adj.bound'))
})

test_that("StratAlphaAdjBSSR certifies the adjusted level over the rectangle", {
  grid <- adj()
  cert <- adj(maximize = 'certified')
  expect_true(all(cert$max.TIE.adj.bound <= 0.1))
  expect_true(all(cert$alpha.adj <= grid$alpha.adj + 1e-15))
  # The level certified is the largest that passes: at the reported level the bound holds,
  # and the type I error rate of the design at that level agrees with the result
  tie <- StratTypeIErrorBSSR(Delta.A = 0.35, N1 = c(6, 6), N2 = c(6, 6), omega = 0.5,
                             r = 1, alpha = cert$alpha.adj[1], tar.power = 0.8,
                             Test = 'CMH', ss.method = 'standard', ss.alpha = 0.1,
                             theta = seq(0.1, 0.9, by = 0.2))
  expect_equal(attr(tie, 'max')$TIE[1], cert$max.TIE.adj[1], tolerance = 1e-12)
  expect_identical(attr(cert, 'nuisance'), 'theta')
})

test_that("StratAlphaAdjBSSR supports the adjustment of the re-estimation and checks input", {
  both <- adj(adjust = 'both', step = 0.005)
  expect_true(all(both$max.TIE.adj <= 0.1))
  k <- (0.1 - both$alpha.adj[1]) / 0.005
  expect_equal(k, round(k), tolerance = 1e-9)
  expect_error(adj(N1 = c(6, 6, 6), N2 = c(6, 6, 6)), 'two strata')
  expect_error(adj(step = 0.2), 'step')
  expect_error(adj(tol = 0), 'tol')
})

test_that("StratAlphaAdjBSSR reports 0 with a warning when the bisection finds no level", {
  # With the whole design observed at the interim analysis and N.max equal to its total,
  # the re-estimation keeps the planned sizes. The largest type I error rate of Cochran's
  # test on the grid is 0.0586 at the level 0.005 and 0.0110 at 0.0025, both above
  # alpha = 0.01, so with tol = 0.3 the bisection ends at the level 0
  expect_warning(
    res <- StratAlphaAdjBSSR(Delta.A = 0.3, N1 = c(6, 6), N2 = c(1, 1), omega = 0.5,
                             r = 6, alpha = 0.01, tar.power = 0.8, Test = 'Cochran',
                             ss.method = 'standard', N.max = 14, maximize = 'grid',
                             theta = seq(0.05, 0.95, by = 0.05), tol = 0.3),
    'reported as 0')
  expect_equal(res$alpha.adj, c(0, 0))
  expect_equal(res$max.TIE.adj, c(0, 0))
  expect_true(all(is.na(res$theta.S1.adj)))
})
