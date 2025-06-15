## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 6,
  fig.align = "center"
)

## ----load_packages, message=FALSE---------------------------------------------
library(bbssr)
library(dplyr)
library(ggplot2)

## ----traditional_power--------------------------------------------------------
# Calculate power for a traditional design
power_result <- BinaryPower(
  p1 = 0.5,          # Response rate in treatment group
  p2 = 0.2,          # Response rate in control group
  N1 = 40,           # Sample size in treatment group
  N2 = 40,           # Sample size in control group
  alpha = 0.025,     # One-sided significance level
  Test = 'Fisher'    # Fisher exact test
)

print(paste("Power:", round(power_result, 3)))

## ----sample_size--------------------------------------------------------------
# Calculate required sample size
sample_size_result <- BinarySampleSize(
  p1 = 0.5,           # Expected response rate in treatment group
  p2 = 0.2,           # Expected response rate in control group
  r = 1,              # Allocation ratio (1:1)
  alpha = 0.025,      # One-sided significance level
  tar.power = 0.8,    # Target power
  Test = 'Boschloo'   # Boschloo exact test (most powerful)
)

print(sample_size_result)

## ----bssr_basic---------------------------------------------------------------
# Basic BSSR calculation
bssr_result <- BinaryPowerBSSR(
  asmd.p1 = 0.45,      # Assumed response rate in treatment group
  asmd.p2 = 0.09,      # Assumed response rate in control group
  p = seq(0.1, 0.9, by = 0.1), # Range of pooled response rates
  Delta.A = 0.36,      # Assumed treatment effect (risk difference)
  Delta.T = 0.36,      # True treatment effect
  N1 = 24,             # Initial sample size in treatment group
  N2 = 24,             # Initial sample size in control group
  omega = 0.5,         # Fraction of data for interim analysis
  r = 1,               # Allocation ratio
  alpha = 0.025,       # Significance level
  tar.power = 0.8,     # Target power
  Test = 'Z-pool',     # Statistical test
  restricted = FALSE,  # Unrestricted design
  weighted = FALSE     # Non-weighted approach
)

# Display first few rows
head(bssr_result)

## ----bssr_comparison, fig.width=11, fig.height=5.5, out.width="100%"----------
# Calculate power of the BSSR with binary endpoint
power.BSSR <- tibble(
  Rule = factor(
    c('Restricted', 'Unrestricted', 'Weighted'), 
    levels = c('Restricted', 'Unrestricted', 'Weighted')
  ),
  restricted = c(TRUE, FALSE, FALSE),
  weighted = c(FALSE, FALSE, TRUE)
) %>% 
  group_by_all() %>% 
  reframe(r = c(1, 2), N1 = c(24, 36), N2 = c(24, 18)) %>% 
  group_by_all() %>% 
  reframe(
    BinaryPowerBSSR(
      asmd.p1 = 0.45, asmd.p2 = 0.09, p = seq(0, 1, by = 0.01),
      Delta.A = 0.36, Delta.T = 0.36, N1, N2, omega = 0.5, r, 
      alpha = 0.025, tar.power = 0.8, Test = 'Z-pool', 
      restricted, weighted
    )
  ) %>% 
  mutate(
    Rule = factor(Rule, levels = c('Restricted', 'Unrestricted', 'Weighted')),
    r_label = paste0('Allocation ratio = ', r, ':1')
  ) %>% 
  rename(Theta = p, Power = power.BSSR)

# Create the improved figure
power.BSSR %>% 
  ggplot(aes(x = Theta, y = Power)) +
  geom_line(aes(color = Rule), linewidth = 1.0) +
  theme_bw() +
  facet_grid(. ~ r_label) + 
  geom_hline(yintercept = 0.8, color = 'gray', linetype = 'longdash', linewidth = 1.0) +
  scale_x_continuous(
    limits = c(0.1, 0.9),
    breaks = seq(0.2, 0.8, by = 0.2),  # Fewer x-axis breaks
    labels = c("0.2", "0.4", "0.6", "0.8")
  ) +
  scale_y_continuous(
    limits = c(0.7, 1.0),
    breaks = seq(0.7, 1.0, by = 0.1),
    labels = c("0.7", "0.8", "0.9", "1.0")
  ) +
  scale_color_manual(
    values = c("Restricted" = "#E31A1C", "Unrestricted" = "#1F78B4", "Weighted" = "#33A02C")
  ) +
  labs(
    x = expression(theta),
    y = "Power",
    title = "Power Comparison: BSSR Design Rules",
    subtitle = "Dashed line indicates target power = 0.8"
  ) +
  theme(
    text = element_text(size = 11),
    plot.title = element_text(size = 13, hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 10, hjust = 0.5, margin = margin(b = 10)),
    panel.spacing = unit(0.8, 'lines'),
    legend.position = 'bottom',
    legend.key.width = unit(1.2, 'cm'),
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.margin = margin(t = 8, r = 0, b = 0, l = 0),
    axis.title.x = element_text(size = 11, margin = margin(t = 8)),
    axis.title.y = element_text(size = 11, margin = margin(r = 8)),
    axis.text = element_text(size = 9),
    strip.text = element_text(size = 10),
    plot.margin = margin(t = 10, r = 10, b = 5, l = 5)
  )

## ----test_comparison----------------------------------------------------------
# Compare different statistical tests
tests <- c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')
test_results <- sapply(tests, function(test) {
  BinaryPower(p1 = 0.5, p2 = 0.2, N1 = 30, N2 = 30, alpha = 0.025, Test = test)
})

test_comparison_df <- data.frame(
  Test = tests,
  Power = round(test_results, 3)
)

print(test_comparison_df)

