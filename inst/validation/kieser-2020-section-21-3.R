# Examines numerically the statement of Kieser (2020, Sect. 21.3, pp. 236-237) that the
# type I error rate of an internal pilot study design is controlled for any blinded
# sample size recalculation rule if the Fisher-Boschloo test is applied in the analysis,
# because the test, as Fisher's exact test on which it is based, conditions on the total
# number of observed events. The statement reports no number, so the quantities that
# bear on it are reported as INFO. Properties that hold by construction, and identities
# between two separate computations of the same quantity, are reported as PASS or FAIL.
#
# Run from the package root after devtools::load_all(); the script uses internal functions
# of the package. The results are written to validation-output/:
#   kieser-2020-section-21-3.csv       every computed quantity with its verdict
#   kieser-2020-section-21-3-rule.csv  the re-estimation rule with the largest type I
#                                      error rate, for each test
#   kieser-2020-section-21-3.md        generated from the two csv files

out.dir <- 'validation-output'
dir.create(out.dir, showWarnings = FALSE)
t.start <- Sys.time()
tab <- list()
add <- function(design, test, item, value, verdict = 'INFO', note = '') {
  tab[[length(tab) + 1L]] <<- data.frame(design = design, test = test, item = item,
                                         value = value, verdict = verdict, note = note,
                                         stringsAsFactors = FALSE)
}
theta.grid <- seq(0.005, 0.995, by = 0.005)

# Designs of Kieser (2020, Example 21.1) as in inst/reproduce/reproduce-published.R:
# one-sided level 0.025, power 0.8, sample size from formula (21.3), unrestricted
# recalculation. Group 1 has the lower rate, so the alternative is 'less'
designs <- list(
  'Delta = 0.30' = list(Delta.A = -0.30, N1 = 39, N2 = 39, n.interim = c(20, 20)),
  'Delta = 0.15' = list(Delta.A = -0.15, N1 = 157, N2 = 157, n.interim = c(79, 79))
)
common <- list(r = 1, alpha = 0.025, tar.power = 0.8, alternative = 'less',
               ss.method = 'standard')
# The refined maximization of the Boschloo p-value is applied to the smaller design only,
# since it is slow for the final sample sizes of the larger one
tests <- list(
  list(label = 'Fisher', Test = 'Fisher', ref.pvalue = FALSE, designs = names(designs)),
  list(label = 'Boschloo', Test = 'Boschloo', ref.pvalue = FALSE,
       designs = names(designs)),
  list(label = 'Boschloo, ref.pvalue', Test = 'Boschloo', ref.pvalue = TRUE,
       designs = 'Delta = 0.30')
)

# Conditional rejection probabilities and the type I error rate of each design --------
for (tt in tests) {
  for (d in tt$designs) {
    args <- c(designs[[d]], common, list(Test = tt$Test, ref.pvalue = tt$ref.pvalue))
    crp <- do.call(BinaryCondRejectBSSR, c(args, list(theta = theta.grid)))
    tie <- do.call(BinaryTypeIErrorBSSR, c(args, list(theta = theta.grid)))
    dif <- max(abs(attr(crp, 'TIE')$TIE - tie$TIE.BSSR))
    add(d, tt$label,
        'difference between the type I error rates from CRP and BinaryTypeIErrorBSSR',
        dif, if (dif < 1e-12) 'PASS' else 'FAIL',
        paste('Largest absolute difference over the grid of theta; the two',
              'computations share only the final sample sizes and the rejection regions'))
    m <- attr(crp, 'max')
    ex <- attr(crp, 'exceed')
    add(d, tt$label, 'number of outcomes (s, s2)', nrow(crp))
    if (tt$Test == 'Fisher') {
      add(d, tt$label, 'largest CRP.total', m$value[2],
          if (m$value[2] <= common$alpha + 1e-12) 'PASS' else 'FAIL',
          paste('Conditional size of Fisher\'s exact test given the total, at most',
                'alpha by construction'))
    } else {
      add(d, tt$label, 'largest CRP.total', m$value[2], 'INFO',
          sprintf('At s = %d, s2 = %d', m$s[2], m$s2[2]))
    }
    add(d, tt$label, 'outcomes with CRP.total above alpha', ex[['CRP.total']])
    add(d, tt$label, 'largest CRP', m$value[1], 'INFO',
        sprintf('At s = %d, s2 = %d', m$s[1], m$s2[1]))
    add(d, tt$label, 'outcomes with CRP above alpha', ex[['CRP']])
    mx <- attr(tie, 'max')
    add(d, tt$label, 'largest type I error rate, re-estimation rule of formula (21.3)',
        mx$TIE[1], 'INFO', sprintf('At theta = %.4f', mx$theta[1]))
    add(d, tt$label, 'largest type I error rate, fixed design', mx$TIE[2], 'INFO',
        sprintf('At theta = %.4f', mx$theta[2]))
  }
}

# Re-estimation rule with the largest type I error rate --------------------------------
# A blinded rule chooses the final sample size from the pooled number s of interim
# responders only. For each value of theta, the rule that maximizes the type I error rate
# picks, for every s, the final size with the largest type I error rate given s. The
# final sizes considered keep the allocation 1:1 and run from the interim size to twice
# the initial size of group 2. The rule chosen at the worst value of theta is then fixed,
# and its type I error rate is recomputed by the binomial summation of
# BinaryTypeIErrorBSSR, which does not use CRP
rules <- list()
d <- 'Delta = 0.30'
for (tt in tests) {
  args <- c(designs[[d]], common, list(Test = tt$Test, ref.pvalue = tt$ref.pvalue))
  n11 <- args$n.interim[1]
  n12 <- args$n.interim[2]
  n1 <- n11 + n12
  cand <- n12:(2 * args$N2)
  # Type I error rate given s, with s in rows, theta in columns and the final size of
  # group 2 in the third dimension
  cond <- vapply(cand, function(N2c) {
    tot <- N2c + max(n11, ceiling(args$r * N2c))
    crp <- do.call(BinaryCondRejectBSSR,
                   c(args, list(N.min = tot, N.max = tot, theta = theta.grid)))
    if (!all(crp$N2 == N2c)) stop('the final size of group 2 is not ', N2c)
    matrix(attr(crp, 'by.s')$TIE.s, nrow = n1 + 1L)
  }, matrix(0, n1 + 1L, length(theta.grid)))
  w <- outer(0:n1, theta.grid, function(s, t) dbinom(s, n1, t))
  worst <- colSums(w * apply(cond, c(1, 2), max))
  j <- which.max(worst)
  rule.N2 <- cand[apply(cond[, j, ], 1, which.max)]
  map <- data.frame(s = 0:n1, N1 = pmax(n11, ceiling(args$r * rule.N2)), N2 = rule.N2)
  attr(map, 'n11') <- n11
  attr(map, 'n12') <- n12
  setup <- bssr_setup(map)
  rr.list <- lapply(seq_along(setup$N1), function(k) {
    get_rr(N1 = setup$N1[k], N2 = setup$N2[k], alpha = args$alpha, Test = args$Test,
           alternative = args$alternative, tsmethod = 'minlike', n.grid = 100,
           bb.gamma = 0, ref.pvalue = args$ref.pvalue, margin = 0)
  })
  f <- function(t) bssr_reject(setup, rr.list, t, t)
  tie.rule <- f(theta.grid)
  dif <- abs(tie.rule[j] - worst[j])
  add(d, tt$label,
      paste('difference between the type I error rates of the worst rule from CRP',
            'and by summation'),
      dif, if (dif < 1e-12) 'PASS' else 'FAIL',
      sprintf('At theta = %.3f, where the rule was chosen', theta.grid[j]))
  add(d, tt$label, 'largest type I error rate over all rules on the grid of theta',
      worst[j], 'INFO',
      sprintf('At theta = %.3f; final sizes of group 2 from %d to %d', theta.grid[j],
              min(cand), max(cand)))
  m <- refine_max(f, theta.grid, tie.rule)
  add(d, tt$label, 'largest type I error rate of the worst rule, refined', m$y, 'INFO',
      sprintf('At theta = %.4f', m$x))
  rules[[length(rules) + 1L]] <- data.frame(design = d, test = tt$label, s = map$s,
                                            N1 = map$N1, N2 = map$N2)
}

# Output ---------------------------------------------------------------------------------
tab <- do.call(rbind, tab)
rules <- do.call(rbind, rules)
utils::write.csv(tab, file.path(out.dir, 'kieser-2020-section-21-3.csv'),
                 row.names = FALSE)
utils::write.csv(rules, file.path(out.dir, 'kieser-2020-section-21-3-rule.csv'),
                 row.names = FALSE)
counts <- table(factor(tab$verdict, levels = c('PASS', 'FAIL', 'INFO')))
md <- c(
  '# Numerical examination of Kieser (2020, Sect. 21.3)',
  '',
  sprintf('Generated by inst/validation/kieser-2020-section-21-3.R on %s (bbssr %s).',
          format(Sys.time(), '%Y-%m-%d %H:%M'),
          as.character(utils::packageVersion('bbssr'))),
  sprintf('Run time: %.1f minutes.',
          as.numeric(difftime(Sys.time(), t.start, units = 'mins'))),
  '',
  paste(names(counts), counts, sep = ': ', collapse = ', '),
  '',
  '| Design | Test | Item | Value | Verdict | Note |',
  '|---|---|---|---|---|---|',
  sprintf('| %s | %s | %s | %s | %s | %s |', tab$design, tab$test, tab$item,
          vapply(tab$value, function(v) format(signif(v, 7)), character(1)),
          tab$verdict, tab$note),
  '',
  '## Re-estimation rule with the largest type I error rate',
  '',
  vapply(split(rules, rules$test), function(r) {
    sprintf('- %s, %s: final size of group 2 for s = 0, ..., %d: %s', r$design[1],
            r$test[1], max(r$s), paste(r$N2, collapse = ', '))
  }, character(1))
)
writeLines(md, file.path(out.dir, 'kieser-2020-section-21-3.md'))
