# Regenerates man/figures/bssr_comparison.png, the design comparison shown in the README.
# Run from the package root with source('tools/make-readme-figure.R').

library(bbssr)
library(ggplot2)

p.grid <- seq(0.19, 0.40, by = 0.03)
args <- list(
  p = p.grid,
  Delta.A = 0.36, Delta.T = 0.36, N1 = 24, N2 = 24, omega = 0.5, r = 1,
  alpha = 0.025, tar.power = 0.8, Test = 'Z-pool'
)
unrestricted <- do.call(BinaryPowerBSSR, c(args, list(restricted = FALSE)))
restricted <- do.call(BinaryPowerBSSR, c(args, list(restricted = TRUE)))

power.df <- data.frame(
  p = rep(unrestricted$p, 3),
  Power = c(unrestricted$power.BSSR, restricted$power.BSSR, unrestricted$power.TRAD),
  Design = factor(
    rep(c('BSSR, unrestricted', 'BSSR, restricted', 'Fixed sample'),
        each = nrow(unrestricted)),
    levels = c('BSSR, unrestricted', 'BSSR, restricted', 'Fixed sample')
  )
)
size.df <- data.frame(
  p = rep(unrestricted$p, 3),
  N = c(unrestricted$E.N, restricted$E.N, rep(48, nrow(unrestricted))),
  Design = power.df$Design
)

palette <- c('BSSR, unrestricted' = 'steelblue',
             'BSSR, restricted' = 'darkorange',
             'Fixed sample' = 'firebrick')

p1 <- ggplot(power.df, aes(x = p, y = Power, colour = Design)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.8, linetype = 'dashed', colour = 'grey40') +
  scale_colour_manual(values = palette) +
  labs(title = 'Power', subtitle = 'Dashed line is the target power of 0.8',
       x = 'True pooled response probability', y = 'Power', colour = NULL) +
  theme_bw() +
  theme(legend.position = 'bottom', plot.title = element_text(face = 'bold'))

p2 <- ggplot(size.df, aes(x = p, y = N, colour = Design)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_colour_manual(values = palette) +
  labs(title = 'Expected total sample size',
       subtitle = 'The fixed design is flat at the planned 48 patients',
       x = 'True pooled response probability', y = 'Expected total N', colour = NULL) +
  theme_bw() +
  theme(legend.position = 'bottom', plot.title = element_text(face = 'bold'))

if (requireNamespace('patchwork', quietly = TRUE)) {
  combined <- patchwork::wrap_plots(p1, p2, ncol = 2, guides = 'collect') &
    theme(legend.position = 'bottom')
  ggsave('man/figures/bssr_comparison.png', combined,
         width = 11, height = 5, dpi = 150)
} else {
  message('patchwork is not installed, saving the power panel only')
  ggsave('man/figures/bssr_comparison.png', p1, width = 7, height = 5, dpi = 150)
}
