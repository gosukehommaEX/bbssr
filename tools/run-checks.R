# Runs every check of the package in one go and writes the results to
# tools/output/run-checks.txt, which is read back during development.
# Open bbssr.Rproj so that the working directory is the package root, then run
#   source('tools/run-checks.R')
# The script needs no object from an earlier session. It takes about five minutes.

out.dir <- file.path('tools', 'output')
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
report <- character(0)
note <- function(...) {
  line <- paste0(...)
  message(line)
  report <<- c(report, line)
}
note('run-checks started ', format(Sys.time(), '%Y-%m-%d %H:%M:%S'))

devtools::document()
devtools::load_all()
note('bbssr version ', as.character(utils::packageVersion('bbssr')))

# Timing on a design grid of five assumed effects and nine interim fractions
timing_grid <- function() {
  grid <- expand.grid(Delta.A = c(0.36, 0.33, 0.30, 0.27, 0.24), omega = 1:9 / 10)
  for (i in seq_len(nrow(grid))) {
    d <- grid$Delta.A[i]
    n0 <- BinarySampleSize(0.4, 0.4 - d, 1, 0.025, 0.9, 'Z-pool')
    BinaryPowerBSSR(p = seq(0.01, 0.99, by = 0.01), Delta.A = d, Delta.T = 0,
                    N1 = n0$N1, N2 = n0$N2, omega = grid$omega[i], r = 1,
                    alpha = 0.025, tar.power = 0.9, Test = 'Z-pool')
  }
}
clear_rr_cache()
t1 <- system.time(timing_grid())[['elapsed']]
t2 <- system.time(timing_grid())[['elapsed']]
note(sprintf('timing: first run %.2f s, second run %.2f s', t1, t2))

# Reproduction of published results, written to reproduce-output/
rep.res <- tryCatch({
  source(file.path('inst', 'reproduce', 'reproduce-published.R'), local = new.env())
  tab <- utils::read.csv(file.path('reproduce-output', 'published-comparison.csv'))
  counts <- table(factor(tab$verdict, levels = c('PASS', 'EXPLAINED', 'FAIL', 'INFO')))
  paste(names(counts), counts, collapse = ', ')
}, error = function(e) paste('ERROR:', conditionMessage(e)))
note('reproduction: ', rep.res)

# Unit tests
test.res <- devtools::test(reporter = 'silent', stop_on_failure = FALSE)
test.df <- as.data.frame(test.res)
note(sprintf('tests: %d expectations, %d failed, %d errors, %d skipped, %d warnings',
             sum(test.df$nb), sum(test.df$failed), sum(test.df$error),
             sum(test.df$skipped), sum(test.df$warning)))
for (t in test.res) {
  for (e in t$results) {
    if (inherits(e, c('expectation_failure', 'expectation_error', 'expectation_warning'))) {
      note('  [', class(e)[1], '] ', t$file, ': ', t$test, ': ',
           gsub('\n', ' | ', conditionMessage(e)))
    }
  }
}

# R CMD check
chk <- devtools::check(error_on = 'never', quiet = TRUE)
note(sprintf('check: %d errors, %d warnings, %d notes', length(chk$errors),
             length(chk$warnings), length(chk$notes)))
for (x in c(chk$errors, chk$warnings, chk$notes)) note('  ', gsub('\n', ' | ', x))

note('run-checks finished ', format(Sys.time(), '%Y-%m-%d %H:%M:%S'))
writeLines(report, file.path(out.dir, 'run-checks.txt'), useBytes = TRUE)
