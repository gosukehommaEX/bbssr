# Computes with bbssr the values plotted in the published figures on matched pairs that the
# vignette 'Validation by reproducing published figures' redraws: Figure 2 of Berger and
# Sidik (2003), Figures 1 and 2 of Fagerland, Lydersen and Laake (2013) and Figures 1 and 2
# of Schepers et al. (2026). The figures themselves are drawn by the functions in
# inst/reproduce/plot-figures.R. The figures of two independent groups are computed by
# inst/reproduce/reproduce-figures.R, which takes about half an hour; this script is kept
# separate so that either can be run alone.
#
# Run from the package root after devtools::load_all(); the conditional power of Schepers
# et al. (2026) uses internal functions of the package. The values are written as CSV files
# to inst/extdata/published-figures/, from which the vignette draws the figures. The time
# of every figure is printed.

out.dir <- file.path('inst', 'extdata', 'published-figures')
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
t.start <- Sys.time()
elapsed <- function(label) {
  message(sprintf('%s: %.1f minutes', label,
                  as.numeric(difftime(Sys.time(), t.start, units = 'mins'))))
}
save_csv <- function(x, name) {
  write.csv(x, file.path(out.dir, name), row.names = FALSE)
}

# Berger and Sidik (2003), Figure 2 ------------------------------------------------------
# Size functions f(p), the rejection probabilities at p01 = p10 = p, of the tests CB (exact
# conditional), M (asymptotic McNemar), Z (exact unconditional) and Z_C (Z with the
# Berger-Boos procedure, beta = 0.0005) for n = 50 and 100 at the one-sided level 0.05,
# with the refinement of the maximization used in reproduce-published.R for Tables 3 to 5
# of the article. The likelihood ratio tests L and L_C of the article are not available in
# bbssr
bs.p <- seq(0, 0.5, by = 0.0025)
bs <- do.call(rbind, lapply(c(50, 100), function(n) {
  size <- function(Test, ...) PairedPower(bs.p, bs.p, n, 0.05, Test, ...)$Power
  data.frame(n = n, p = bs.p, CB = size('McNemar-exact'), M = size('McNemar'),
             Z = size('McNemar-uncond', ref.pvalue = TRUE),
             Z.C = size('McNemar-uncond', ref.pvalue = TRUE, bb.gamma = 5e-4))
}))
save_csv(bs, 'berger-sidik-2003-figure2.csv')
elapsed('Berger and Sidik (2003), Figure 2')

# Fagerland, Lydersen and Laake (2013) ---------------------------------------------------
# The article gives the cell probabilities by the marginal success probabilities p1+ and
# p+1 and the odds ratio theta = p11 p22 / (p12 p21). p11 is the root of
# p11 p22 = theta p12 p21 that lies within the bounds given by the margins, which the
# article does not print; with it the type I error rates of Table 4 of the article are
# reproduced in reproduce-published.R. The discordant cells then have the probabilities
# p12 = p1+ - p11 and p21 = p+1 - p11. The closed form below is the root of the quadratic
# (1 - theta) x^2 + (1 + (theta - 1) (p1 + p2)) x - theta p1 p2 = 0 written without
# cancellation, for theta >= 1
fll_p11 <- function(p1, p2, theta) {
  b <- 1 + (theta - 1) * (p1 + p2)
  2 * theta * p1 * p2 / (b + sqrt(b^2 + 4 * theta * (1 - theta) * p1 * p2))
}
# The five tests of the article at the two-sided level 0.05. The two-sided exact
# unconditional p-value of the article is twice the one-sided p-value, which is the p-value
# of bbssr ordered by the absolute value of the statistic, since the null distribution of
# the outcomes is symmetric
fll.tests <- c(asymptotic = 'McNemar', cc = 'McNemar-cc', exact = 'McNemar-exact',
               midp = 'McNemar-midP', uncond = 'McNemar-uncond')
fll_rates <- function(p10, p01, N) {
  vapply(fll.tests, function(tst) {
    PairedPower(p10, p01, N, 0.05, tst, 'two.sided',
                ref.pvalue = tst == 'McNemar-uncond')$Power
  }, numeric(length(p10)))
}

# Figure 1: type I error rates against p1+ = p+1 for (N, theta) = (20, 2), (40, 1),
# (70, 10) and (100, 3)
fll.p <- seq(0, 1, by = 0.005)
fll1 <- do.call(rbind, lapply(list(c(20, 2), c(40, 1), c(70, 10), c(100, 3)), function(s) {
  q <- pmax(fll.p - fll_p11(fll.p, fll.p, s[2]), 0)
  data.frame(N = s[1], theta = s[2], p = fll.p, fll_rates(q, q, s[1]))
}))
save_csv(fll1, 'fagerland-2013-figure1.csv')
elapsed('Fagerland, Lydersen and Laake (2013), Figure 1')

# Figure 2: power against the number of matched pairs for p1+ = 0.1, p+1 = 0.35 and
# theta = 2. The two-sided tests do not depend on which discordant cell is p10
fll.p11 <- fll_p11(0.1, 0.35, 2)
fll2 <- data.frame(N = 1:100, t(vapply(1:100, function(N) {
  fll_rates(0.1 - fll.p11, 0.35 - fll.p11, N)
}, numeric(length(fll.tests)))))
save_csv(fll2, 'fagerland-2013-figure2.csv')
elapsed('Fagerland, Lydersen and Laake (2013), Figure 2')

# Schepers et al. (2026), Figures 1 and 2 -----------------------------------------------
# Two-sided exact conditional McNemar test at level 0.05 with n(1) = 38 interim pairs, of
# which nd(1) = 0, ..., 38 are discordant, and p* = 0.50, 0.51, ..., 0.90, the probability
# that a discordant pair favours treatment 1. The conditional power (4) of the article
# needs the probability psi of a discordant pair in the second stage. The figures agree
# with the planning value psi = 0.2 of the article and not with the interim estimate
# nd(1) / n(1), which the re-estimation rule of Section 2.4 uses; the conditional power
# with the estimate is stored as well. The article does not state the largest
# second-stage size of Figure 2: its unattainable cells agree with 70, and the smallest
# second-stage size is stored for sizes up to 100 so that the effect of the bound can be
# seen
sc.ps <- round(seq(0.5, 0.9, by = 0.01), 2)
sc.nd <- 0:38
# Cells ordered by p* first, so that paired_cp() computes the rejection probabilities
# given the number of discordant pairs once for each p*
sc.cells <- data.frame(nd = rep(sc.nd, times = length(sc.ps)),
                       p.star = rep(sc.ps, each = length(sc.nd)))
sc_cp <- function(n2, psi) {
  rr <- bbssr:::get_paired_rr(38 + n2, 0.05, 'McNemar-exact', 'two.sided', 100, 0, FALSE, 0)
  bbssr:::paired_cp(rr, 38L, sc.cells$nd, psi, sc.cells$p.star)
}
# Figure 1: conditional power with the fixed second-stage size n(2) = 38
sc1 <- data.frame(sc.cells, CP = sc_cp(38, rep(0.2, nrow(sc.cells))),
                  CP.estimate = sc_cp(38, sc.cells$nd / 38))
save_csv(sc1, 'schepers-2026-figure1.csv')
elapsed('Schepers et al. (2026), Figure 1')

# Figure 2: smallest second-stage size n(2) = 0, ..., 100 with a conditional power of at
# least 0.8, NA if there is none
sc2 <- data.frame(sc.cells, n2.min = NA_integer_)
for (n2 in 0:100) {
  hit <- is.na(sc2$n2.min) & sc_cp(n2, rep(0.2, nrow(sc.cells))) >= 0.8
  sc2$n2.min[hit] <- n2
}
save_csv(sc2, 'schepers-2026-figure2.csv')
elapsed('Schepers et al. (2026), Figure 2')
