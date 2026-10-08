# Design of a small re-estimation trial with the rejection regions of its final sample
# sizes and the index of the final sample size reached from each pooled interim count
crp_design <- function(Test, alternative, r, n.interim, Delta.A) {
  map <- bssr_map(Delta.A = Delta.A, N1 = 10, N2 = 10, omega = NULL,
                  n.interim = n.interim, r = r, alpha = 0.025, tar.power = 0.8,
                  Test = Test, restricted = FALSE, alternative = alternative,
                  tsmethod = 'minlike', n.grid = 100L, bb.gamma = 0, effect = 'RD',
                  ss.method = 'standard', ss.Test = Test, ss.alpha = 0.025,
                  rounding = 'group', N.min = NULL, N.max = NULL, ref.pvalue = FALSE,
                  margin = 0, margin.scale = 'RD')
  st <- bssr_setup(map)
  rr.list <- lapply(seq_along(st$N1), function(k) {
    get_rr(N1 = st$N1[k], N2 = st$N2[k], alpha = 0.025, Test = Test,
           alternative = alternative, tsmethod = 'minlike', n.grid = 100, bb.gamma = 0,
           ref.pvalue = FALSE, margin = 0, margin.scale = 'RD')
  })
  id <- match(paste(map$N1, map$N2), paste(st$N1, st$N2)) - 1L
  list(st = st, rr.list = rr.list, id = as.integer(id))
}

test_that("bssr_cond_reject agrees with a direct summation", {
  cases <- list(
    list(Test = 'Fisher', alternative = 'greater', r = 1, n.interim = c(4, 3),
         Delta.A = 0.4),
    list(Test = 'Chisq', alternative = 'two.sided', r = 2, n.interim = c(4, 2),
         Delta.A = 0.4)
  )
  for (cs in cases) {
    d <- do.call(crp_design, cs)
    got <- bssr_cond_reject(d$rr.list, d$id, d$st$n21, d$st$n22, d$st$n11, d$st$n12)
    want <- bssr_cond_reject_ref(d$rr.list, d$id, d$st$n21, d$st$n22, d$st$n11,
                                 d$st$n12)
    expect_gt(length(got$crp), 0)
    expect_equal(got$s, want$s)
    expect_equal(got$s2, want$s2)
    expect_equal(got$crp, want$crp, tolerance = 1e-12)
    expect_equal(got$crp.total, want$crp.total, tolerance = 1e-12)
  }
})

test_that("bssr_cond_reject rejects inconsistent input", {
  d <- crp_design(Test = 'Fisher', alternative = 'greater', r = 1, n.interim = c(4, 3),
                  Delta.A = 0.4)
  st <- d$st
  expect_error(bssr_cond_reject(d$rr.list, d$id[-1], st$n21, st$n22, st$n11, st$n12),
               'one element per pooled number')
  bad <- d$id
  bad[1] <- length(d$rr.list)
  expect_error(bssr_cond_reject(d$rr.list, bad, st$n21, st$n22, st$n11, st$n12),
               'outside rr_list')
  expect_error(bssr_cond_reject(d$rr.list, d$id, st$n21 + 1L, st$n22, st$n11, st$n12),
               'dimension')
  expect_error(bssr_cond_reject(d$rr.list, d$id, st$n21[-1], st$n22, st$n11, st$n12),
               'one element per rejection region')
})
