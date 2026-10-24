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

test_that("StratPower of the unstratified tests agrees with the reference values", {
  both <- function(Test, alternative) c(pw2(Test, alternative), pw3(Test, alternative))
  expect_equal(both('Unstrat-Chisq', 'greater'), c(0.335288560357361, 0.177197685152934),
               tolerance = 1e-10)
  expect_equal(both('Unstrat-Chisq', 'less'), c(0.00202571746993205, 0.025213285248676),
               tolerance = 1e-10)
  expect_equal(both('Unstrat-Chisq', 'two.sided'), c(0.210566268355524, 0.15294743714339),
               tolerance = 1e-10)
  expect_equal(both('Unstrat-Chisq-cc', 'greater'),
               c(0.208522344468668, 0.0970165390611851),
               tolerance = 1e-10)
  expect_equal(both('Unstrat-Chisq-cc', 'less'),
               c(0.000644444554309201, 0.0101001563835337),
               tolerance = 1e-10)
  expect_equal(both('Unstrat-Chisq-cc', 'two.sided'),
               c(0.115119583787729, 0.0535554775819315),
               tolerance = 1e-10)
  expect_equal(both('Unstrat-Fisher', 'greater'), c(0.208522344468668, 0.0970165390611851),
               tolerance = 1e-10)
  expect_equal(both('Unstrat-Fisher', 'less'), c(0.000644444554309201, 0.0101001563835337),
               tolerance = 1e-10)
  expect_equal(both('Unstrat-Fisher', 'two.sided'),
               c(0.115119583787729, 0.0828363802829164),
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

test_that("StratPower reproduces the design example and Table 4 of Jung (2014)", {
  # Design example (p. 139): three strata with the control probabilities
  # (0.9, 0.75, 0.6), the odds ratios (1, 30, 30), 10, 10 and 11 patients per group and
  # the one-sided level 0.1: marginal type I error rate 0.0565 and power 0.8105
  q <- c(0.9, 0.75, 0.6)
  p <- q * c(1, 30, 30) / (1 - q + q * c(1, 30, 30))
  m <- c(10, 10, 11)
  expect_equal(round(StratPower(q, q, m, m, 0.1, 'Strat-Fisher')$Power, 4), 0.0565)
  expect_equal(round(StratPower(p, q, m, m, 0.1, 'Strat-Fisher')$Power, 4), 0.8105)
  # Table 4, both allocations fixed: the smallest N from the size of the Mantel-Haenszel
  # test upwards with power 0.9 at the one-sided level 0.05, with the stratum fractions
  # (0.25, 0.75), the fractions b of group 1 and the odds ratios (5, 10)
  jung_n <- function(b, N) {
    q <- c(0.1, 0.3)
    p <- q * c(5, 10) / (1 - q + q * c(5, 10))
    repeat {
      n1 <- floor(N * 0.25)
      m <- floor(N * c(0.25, 0.75) * b)
      if (StratPower(p, q, m, c(n1, N - n1) - m, 0.05, 'Strat-Fisher')$Power >= 0.9) {
        return(N)
      }
      N <- N + 1
    }
  }
  expect_equal(c(jung_n(c(0.25, 0.25), 46), jung_n(c(0.5, 0.5), 36)), c(59, 43))
})

test_that("StratPower with a random allocation agrees with the reference values", {
  # Reference values from tools/reference/reference_strat.py (alloc_power), which sums the
  # power of the fixed design of every configuration of stratum and group sizes: N = 16,
  # equal fractions, r = 1, blocks of 4, one-sided level 0.05; the first scenario has an
  # effect and the second none
  pw <- function(alloc, Test) {
    StratPower(rbind(c(0.5, 0.8), c(0.3, 0.6)), rbind(c(0.2, 0.4), c(0.3, 0.6)),
               alpha = 0.05, Test = Test, alloc = alloc, N = 16, frac = c(0.5, 0.5))$Power
  }
  expect_equal(pw('stratum', 'CMH'), c(0.411444606689659, 0.0473474558748814),
               tolerance = 1e-10)
  expect_equal(pw('stratum', 'Strat-Fisher'), c(0.223592332438072, 0.0142843045203167),
               tolerance = 1e-10)
  expect_equal(pw('stratum', 'Unstrat-Chisq'), c(0.383824090759094, 0.0394443754222324),
               tolerance = 1e-10)
  expect_equal(pw('group', 'CMH'), c(0.396656774694748, 0.0501045228602838),
               tolerance = 1e-10)
  expect_equal(pw('group', 'Strat-Fisher'), c(0.199885429075408, 0.0131489604395502),
               tolerance = 1e-10)
  expect_equal(pw('group', 'Unstrat-Chisq'), c(0.405920505935753, 0.0544615541355178),
               tolerance = 1e-10)
  expect_equal(pw('both', 'CMH'), c(0.396610240208928, 0.0500683414440941),
               tolerance = 1e-10)
  expect_equal(pw('both', 'Strat-Fisher'), c(0.199358926866081, 0.0131096992312949),
               tolerance = 1e-10)
  expect_equal(pw('both', 'Unstrat-Chisq'), c(0.407803164480989, 0.0545926491870062),
               tolerance = 1e-10)
  # Fractions (0.4, 0.6) and the allocation ratio 3: stratum sizes (6, 11) with 'group',
  # and blocks of 4 holding 3 patients of group 1 with 'stratum'
  pu <- function(alloc) {
    StratPower(c(0.5, 0.8), c(0.2, 0.4), alpha = 0.05, Test = 'Strat-Fisher', alloc = alloc,
               N = 17, frac = c(0.4, 0.6), r = 3)
  }
  expect_equal(c(pu('stratum')$Power, pu('group')$Power),
               c(0.190973505647753, 0.16684528129529),
               tolerance = 1e-10)
  expect_identical(attr(pu('group'), 'n.strata'), c(6L, 11L))
  expect_null(attr(pu('stratum'), 'n.strata'))
  expect_identical(attr(pu('stratum'), 'block'), 4L)
})

test_that("StratPower with fixed stratum sizes reproduces Table 4 of Jung (2014)", {
  # Table 4, column n fixed (Section 3.2): stratum sizes floor(N / 2) and the rest, group 1
  # sizes binomial with b = (0.5, 0.5), (q1, q2) = (0.1, 0.3), odds ratios (5, 10),
  # stratified exact test at the one-sided level 0.05. The table gives N = 56, the first
  # size with power 0.9; the powers 0.8998 at N = 55 and 0.9032 at N = 56 were also
  # computed independently
  q <- c(0.1, 0.3)
  p <- q * c(5, 10) / (1 - q + q * c(5, 10))
  pw <- function(N) {
    StratPower(p, q, alpha = 0.05, Test = 'Strat-Fisher', alloc = 'group', N = N,
               frac = c(0.5, 0.5))$Power
  }
  expect_equal(round(c(pw(55), pw(56)), 4), c(0.8998, 0.9032))
})

test_that("StratPower with a random allocation handles empty strata and groups", {
  # With two patients every configuration leaves a stratum or a group empty, and no test
  # can reject at the level 0.05
  for (alloc in c('stratum', 'group', 'both')) {
    for (Test in c('CMH', 'Strat-Fisher', 'Unstrat-Chisq')) {
      res <- StratPower(c(0.5, 0.8), c(0.2, 0.4), alpha = 0.05, Test = Test, alloc = alloc,
                        N = 2, frac = c(0.5, 0.5), block = 2)
      expect_identical(res$Power, 0, info = paste(alloc, Test))
    }
  }
})

test_that("StratPower checks the arguments of the allocation", {
  expect_error(StratPower(c(0.5, 0.7), c(0.3, 0.5), alpha = 0.05, Test = 'CMH'),
               'N1 and N2 must be supplied')
  expect_error(StratPower(c(0.5, 0.7), c(0.3, 0.5), c(8, 6), c(7, 9), 0.05, 'CMH', N = 20),
               'random allocation')
  expect_error(StratPower(c(0.5, 0.7), c(0.3, 0.5), c(8, 6), c(7, 9), 0.05, 'CMH',
                          alloc = 'both', N = 20, frac = c(1, 1)), 'not by N1')
  expect_error(StratPower(c(0.5, 0.7), c(0.3, 0.5), alpha = 0.05, Test = 'CMH',
                          alloc = 'both', frac = c(1, 1)), 'N must be supplied')
  expect_error(StratPower(c(0.5, 0.7, 0.4), c(0.3, 0.5, 0.2), alpha = 0.05, Test = 'CMH',
                          alloc = 'both', N = 20, frac = c(1, 1, 1)), 'two strata')
})
