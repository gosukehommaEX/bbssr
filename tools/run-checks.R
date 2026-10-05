# Runs the checks of the package in one go and writes the results to
# tools/output/run-checks.txt, which is read back during development.
# Open bbssr.Rproj so that the working directory is the package root, then run
#   source('tools/run-checks.R')
# for the quick run (documentation, timing of the load_all() build and unit tests),
#   reproduce <- TRUE; source('tools/run-checks.R')
# to add the reproduction of published results and the validation scripts, or
#   full.check <- TRUE; source('tools/run-checks.R')
# to run everything, including the timing of an optimized installed build and R CMD
# check. The flags are removed when the script starts, so the next run is quick again.
# The script needs no other object from an earlier session.

flag <- function(name) {
  on <- exists(name, envir = globalenv()) && isTRUE(get(name, envir = globalenv()))
  if (exists(name, envir = globalenv())) rm(list = name, envir = globalenv())
  on
}
run.full <- flag('full.check')
run.reproduce <- flag('reproduce') || run.full
out.dir <- file.path('tools', 'output')
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
report <- character(0)
# Each line is written to the report at once, so an interrupted run keeps what it found
note <- function(...) {
  line <- paste0(...)
  message(line)
  report <<- c(report, line)
  writeLines(report, file.path(out.dir, 'run-checks.txt'), useBytes = TRUE)
}
note('run-checks started ', format(Sys.time(), '%Y-%m-%d %H:%M:%S'),
     if (run.full) ' (full run)' else if (run.reproduce) ' (quick run with reproduction)'
     else ' (quick run)')

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
note(sprintf('timing, load_all() build without optimization: first run %.2f s, second run %.2f s',
             t1, t2))
# The same timing with an optimized build, installed into a temporary library and run in
# a separate R process, which is how the package is used after installation. The
# installation runs R CMD INSTALL in its own process, because install.packages() on
# Windows refuses to install a package that is loaded in the current session
opt <- if (!run.full) 'skipped in the quick run' else tryCatch({
  lib <- file.path(tempdir(), 'bbssr-timing-lib')
  dir.create(lib, showWarnings = FALSE)
  tarball <- pkgbuild::build(dest_path = tempdir(), vignettes = FALSE, quiet = TRUE)
  inst <- callr::rcmd('INSTALL', c(paste0('--library=', lib), tarball), fail_on_status = FALSE)
  if (inst$status != 0) stop('R CMD INSTALL failed: ', inst$stderr)
  callr::r(function(lib, f) {
    library(bbssr, lib.loc = lib)
    environment(f) <- globalenv()
    c(system.time(f())[['elapsed']], system.time(f())[['elapsed']])
  }, args = list(lib = lib, f = timing_grid))
}, error = function(e) cli::ansi_strip(conditionMessage(e)))
if (is.numeric(opt)) {
  note(sprintf('timing, optimized installed build: first run %.2f s, second run %.2f s',
               opt[1], opt[2]))
} else if (run.full) {
  note('timing, optimized installed build: ERROR ', opt)
} else {
  note('timing, optimized installed build: ', opt)
}

# Reproduction of published results, written to reproduce-output/
rep.res <- if (!run.reproduce) 'skipped' else tryCatch({
  source(file.path('inst', 'reproduce', 'reproduce-published.R'), local = new.env())
  tab <- utils::read.csv(file.path('reproduce-output', 'published-comparison.csv'))
  counts <- table(factor(tab$verdict, levels = c('PASS', 'EXPLAINED', 'FAIL', 'INFO')))
  paste(names(counts), counts, collapse = ', ')
}, error = function(e) paste('ERROR:', cli::ansi_strip(conditionMessage(e))))
note('reproduction: ', rep.res)

# Numerical examination of Kieser (2020, Sect. 21.3), written to validation-output/
val.res <- if (!run.reproduce) 'skipped' else tryCatch({
  source(file.path('inst', 'validation', 'kieser-2020-section-21-3.R'), local = new.env())
  tab <- utils::read.csv(file.path('validation-output', 'kieser-2020-section-21-3.csv'))
  counts <- table(factor(tab$verdict, levels = c('PASS', 'FAIL', 'INFO')))
  paste(names(counts), counts, collapse = ', ')
}, error = function(e) paste('ERROR:', cli::ansi_strip(conditionMessage(e))))
note('validation: ', val.res)

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
           gsub('\n', ' | ', cli::ansi_strip(conditionMessage(e))))
    }
  }
}

# R CMD check
if (run.full) {
  chk <- devtools::check(error_on = 'never', quiet = TRUE)
  note(sprintf('check: %d errors, %d warnings, %d notes', length(chk$errors),
               length(chk$warnings), length(chk$notes)))
  for (x in c(chk$errors, chk$warnings, chk$notes)) note('  ', gsub('\n', ' | ', x))
} else {
  note('check: skipped in the quick run')
}

note('run-checks finished ', format(Sys.time(), '%Y-%m-%d %H:%M:%S'))
writeLines(report, file.path(out.dir, 'run-checks.txt'), useBytes = TRUE)
