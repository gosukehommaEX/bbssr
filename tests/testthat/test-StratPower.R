pw2 <- function(Test, alternative) {
  StratPower(c(0.6, 0.45), c(0.35, 0.3), c(8, 6), c(7, 9), 0.05, Test, alternative,
             'central')$Power
}
pw3 <- function(Test, alternative) {
  StratPower(c(0.7, 0.3, 0.55), c(0.4, 0.25, 0.5), c(4, 5, 3), c(5, 3, 4), 0.1, Test,
             alternative, 'central')$Power
}

test_that("StratPower agrees with the reference values with two and three strata", {
  both <- function(Test, alternative) c(pw2(Test, alternative), pw3(Test, alternative))
  expect_equal(both('CMH', 'greater'), c(0.303344073443283, 0.288160817883127),
               tolerance = 1e-10)
  expect_equal(both('CMH', 'less'), c(0.00281539830312937, 0.0231946703751953),
               tolerance = 1e-10)
  expect_equal(both('CMH', 'two.sided'), c(0.21361470669462, 0.183714391259765),
               tolerance = 1e-10)
  expect_equal(both('CMH-cc', 'greater'), c(0.175986004477007, 0.157304458861354),
               tolerance = 1e-10)
  expect_equal(both('CMH-cc', 'less'), c(0.000750915604308426, 0.0077223886922672),
               tolerance = 1e-10)
  expect_equal(both('CMH-cc', 'two.sided'), c(0.116205630058535, 0.0810268845923619),
               tolerance = 1e-10)
  expect_equal(both('Cochran', 'greater'), c(0.340718272696192, 0.322681717131867),
               tolerance = 1e-10)
  expect_equal(both('Cochran', 'less'), c(0.00344534812589269, 0.0264725126160019),
               tolerance = 1e-10)
  expect_equal(both('Cochran', 'two.sided'), c(0.215932581835288, 0.21722570075209),
               tolerance = 1e-10)
  expect_equal(both('Cochran-cc', 'greater'), c(0.211447975910345, 0.177704777618129),
               tolerance = 1e-10)
  expect_equal(both('Cochran-cc', 'less'), c(0.000912709420274477, 0.00836511632673804),
               tolerance = 1e-10)
  expect_equal(both('Cochran-cc', 'two.sided'), c(0.1204093750041, 0.104361445463374),
               tolerance = 1e-10)
  expect_equal(both('Strat-Fisher', 'greater'), c(0.199384613187037, 0.161792083831048),
               tolerance = 1e-10)
  expect_equal(both('Strat-Fisher', 'less'), c(0.000912316604662701, 0.00781548628788537),
               tolerance = 1e-10)
  expect_equal(both('Strat-Fisher', 'two.sided'), c(0.119332788853109, 0.0925061420069531),
               tolerance = 1e-10)
})

test_that("StratPower evaluates several scenarios and recycles a single one", {
  p1 <- rbind(c(0.6, 0.45), c(0.4, 0.4))
  res <- StratPower(p1, c(0.35, 0.3), c(8, 6), c(7, 9), 0.05, 'CMH')
  expect_s3_class(res, 'bbssr_spower')
  expect_named(res, c('p1.S1', 'p1.S2', 'p2.S1', 'p2.S2', 'Power'))
  expect_equal(res$Power[1], pw2('CMH', 'greater'))
  expect_equal(res$p2.S2, c(0.3, 0.3))
  expect_identical(attr(res, 'N1'), c(8L, 6L))
  expect_error(StratPower(p1, rbind(c(0.3, 0.3), c(0.3, 0.3), c(0.3, 0.3)), c(8, 6),
                          c(7, 9), 0.05, 'CMH'), 'same number of scenarios')
  expect_error(StratPower(c(0.6, 0.45, 0.2), c(0.35, 0.3), c(8, 6), c(7, 9), 0.05, 'CMH'),
               'one element per stratum')
})

test_that("StratPower of the stratified exact test does not exceed the level under the null", {
  res <- StratPower(rbind(c(0.2, 0.6), c(0.5, 0.5), c(0.9, 0.1)),
                    rbind(c(0.2, 0.6), c(0.5, 0.5), c(0.9, 0.1)), c(10, 7), c(9, 12), 0.05,
                    'Strat-Fisher', 'two.sided', 'minlike')
  expect_true(all(res$Power <= 0.05))
})
