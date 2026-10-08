# Runs the checks of the package in one go and writes the results to
# tools/output/run-checks.txt, which is read back during development.
# Open bbssr.Rproj so that the working directory is the package root, then run
#   source('tools/run-checks.R')
# for the quick run (documentation, timing of the load_all() build and unit tests),
#   reproduce <- TRUE; source('tools/run-checks.R')
# to add the reproduction of published results and the scripts in tools/validation/,
#   docs <- TRUE; source('tools/run-checks.R')
# to add README.md from README.Rmd, every vignette rendered in a separate R process with
# its time, the spelling check and the pkgdown site in docs/, or
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
run.docs <- flag('docs') || run.full
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
     else if (run.docs) ' (quick run with documentation)' else ' (quick run)')

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

# Every script in tools/validation/, if the folder is present, run from the package root
val.files <- sort(list.files(file.path('tools', 'validation'), pattern = '[.]R$',
                             full.names = TRUE))
if (!run.reproduce) {
  note('validation: skipped')
} else if (length(val.files) == 0) {
  note('validation: no script present')
} else {
  for (f in val.files) {
    val.res <- tryCatch({
      sprintf('done in %.1f s',
              system.time(source(f, local = new.env()))[['elapsed']])
    }, error = function(e) paste('ERROR:', cli::ansi_strip(conditionMessage(e))))
    note('validation, ', basename(f), ': ', val.res)
  }
}

# Unit tests
test.res <- devtools::test(reporter = 'silent', stop_on_failure = FALSE)
test.df <- as.data.frame(test.res)
note(sprintf('tests: %d expectations, %d failed, %d errors, %d skipped, %d warnings',
             sum(test.df$nb), sum(test.df$failed), sum(test.df$error),
             sum(test.df$skipped), sum(test.df$warning)))
# Time of every test_that() block, slowest first, since each block should finish within
# 10 seconds
test.times <- test.df[order(-test.df$real), c('file', 'test', 'real')]
utils::write.csv(test.times, file.path(out.dir, 'test-times.csv'), row.names = FALSE)
note(sprintf('tests: slowest block %.1f s (%s: %s), %d block(s) over 10 s, listed in %s',
             test.times$real[1], test.times$file[1], test.times$test[1],
             sum(test.times$real > 10), file.path(out.dir, 'test-times.csv')))
for (t in test.res) {
  for (e in t$results) {
    if (inherits(e, c('expectation_failure', 'expectation_error', 'expectation_warning'))) {
      note('  [', class(e)[1], '] ', t$file, ': ', t$test, ': ',
           gsub('\n', ' | ', cli::ansi_strip(conditionMessage(e))))
    }
  }
}

# Documentation: README.md from README.Rmd, every vignette rendered in a separate R
# process from the installed package, the spelling check and the pkgdown site in docs/.
# The package is installed into the user library first, since the path of the package
# contains characters that pkgdown cannot install from, and the remaining steps are
# skipped when the installation fails. Warnings raised in a chunk are written into the
# rendered output by knitr, so they are counted there
run.docs.steps <- run.docs
if (run.docs) {
  inst.res <- tryCatch({
    # reload = FALSE keeps the installed DLL out of this session. load_all() does not
    # unload a DLL loaded from the library, and Windows refuses to overwrite a DLL in use,
    # so a reloaded package makes the installation fail in a later run of this session
    devtools::install(upgrade = FALSE, quiet = TRUE, reload = FALSE)
    paste('bbssr', callr::r(function() as.character(utils::packageVersion('bbssr'))),
          'installed')
  }, error = function(e) {
    # The output of R CMD INSTALL is kept in the error object, not in its message
    out <- character(0)
    x <- e
    while (inherits(x, 'condition')) {
      for (f in c('stdout', 'stderr')) {
        if (is.character(x[[f]])) out <- c(out, strsplit(x[[f]], '\r?\n')[[1]])
      }
      x <- x$parent
    }
    out <- utils::tail(cli::ansi_strip(out[nzchar(trimws(out))]), 10)
    paste(c(paste('ERROR:', cli::ansi_strip(conditionMessage(e))), out), collapse = ' | ')
  })
  note('docs, installation: ', inst.res)
  run.docs.steps <- !startsWith(inst.res, 'ERROR')
}
if (run.docs.steps) {
  readme.res <- tryCatch({
    callr::r(function(f) {
      rmarkdown::render(f, output_options = list(html_preview = FALSE), quiet = TRUE,
                        envir = new.env())
    }, args = list(f = normalizePath('README.Rmd')))
    sprintf('README.md written, %d warnings in the output',
            sum(grepl('^#> Warning', readLines('README.md', warn = FALSE))))
  }, error = function(e) paste('ERROR:', cli::ansi_strip(conditionMessage(e))))
  note('docs, README: ', readme.res)
  for (f in sort(list.files('vignettes', pattern = '[.]Rmd$', full.names = TRUE))) {
    vig.res <- tryCatch({
      res <- callr::r(function(f, out) {
        t <- system.time(html <- rmarkdown::render(f, output_dir = out,
                                                   intermediates_dir = out,
                                                   quiet = TRUE, envir = new.env()))
        list(time = t[['elapsed']],
             warnings = sum(grepl('#&gt; Warning', readLines(html, warn = FALSE))))
      }, args = list(f = normalizePath(f), out = tempfile('vignette-')))
      sprintf('%.1f s, %d warnings in the output', res$time, res$warnings)
    }, error = function(e) paste('ERROR:', cli::ansi_strip(conditionMessage(e))))
    note('docs, vignette ', basename(f), ': ', gsub('\n', ' | ', vig.res))
  }
  spell.res <- tryCatch({
    sp <- spelling::spell_check_package()
    words <- if (nrow(sp) == 0) character(0) else {
      vapply(seq_len(nrow(sp)), function(i) {
        paste0(sp$word[i], ': ', paste(sp$found[[i]], collapse = ', '))
      }, character(1))
    }
    writeLines(words, file.path(out.dir, 'spelling.txt'), useBytes = TRUE)
    sprintf('%d words not in the dictionary or inst/WORDLIST, listed in %s',
            nrow(sp), file.path(out.dir, 'spelling.txt'))
  }, error = function(e) paste('ERROR:', cli::ansi_strip(conditionMessage(e))))
  note('docs, spelling: ', spell.res)
  site.res <- tryCatch({
    pkgdown::build_site(install = FALSE, preview = FALSE)
    'site written to docs/'
  }, error = function(e) paste('ERROR:', cli::ansi_strip(conditionMessage(e))))
  note('docs, pkgdown: ', site.res)
} else if (!run.docs) {
  note('docs: skipped')
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
