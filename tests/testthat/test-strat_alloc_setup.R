setup_design <- function(N1) {
  strat_design(N1, N1, 1, NULL, 0.5, NULL, 0.05, 'CMH', 'greater', 'minlike', 'standard',
               1e5)
}
re_args <- function(d) {
  function(level) {
    list(Delta.A = 0.3, effect = 'RD', frac = d$frac, r = d$r, ss.method = 'standard',
         ss.Test = 'CMH', ss.alpha = level, tar.power = 0.8, N1.plan = d$a$N1,
         N2.plan = d$a$N2, restricted = FALSE, N.min = NULL, N.max = NULL,
         search = 'crossing', search.limit = c(2, 50), ss.strata = 'stratified')
  }
}

test_that("strat_alloc_setup keeps the draws of the sequences across the levels", {
  d <- setup_design(c(6, 6))
  set.seed(10)
  st <- strat_alloc_setup(d, 'both', 4, 5, re_args(d), 'observed', FALSE, NULL, 60)
  expect_equal(st$n, 12)
  expect_equal(st$N.plan, 24)
  a <- st$layout(0.05)
  b <- st$layout(0.2)
  expect_length(a$inputs, 5)
  for (i in 1:5) {
    k <- min(nrow(a$inputs[[i]]$counts), nrow(b$inputs[[i]]$counts))
    expect_identical(a$inputs[[i]]$counts[seq_len(k), ], b$inputs[[i]]$counts[seq_len(k), ])
  }
  # A larger level of the re-estimation gives final totals that are not larger
  Na <- unlist(lapply(a$inputs, `[[`, 'Nfin'))
  Nb <- unlist(lapply(b$inputs, `[[`, 'Nfin'))
  expect_true(all(Nb <= Na) && any(Nb < Na))
  expect_true(all(Na >= 12 & Na <= 60))
})

test_that("strat_alloc_setup checks the settings", {
  d <- setup_design(c(6, 6))
  expect_error(strat_alloc_setup(d, 'both', 4, 1, re_args(d), 'observed', FALSE, NULL,
                                 NULL), 'nseq')
  expect_error(strat_alloc_setup(d, 'both', 4, 5, re_args(d), 'observed', FALSE, NULL, 10),
               'N.max')
  d3 <- setup_design(c(6, 6, 6))
  expect_error(strat_alloc_setup(d3, 'both', 4, 5, re_args(d3), 'observed', FALSE, NULL,
                                 NULL), 'two strata')
})
