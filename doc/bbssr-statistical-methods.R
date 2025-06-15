## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 6,
  fig.align = "center"
)

## ----load_packages, message=FALSE, warning=FALSE------------------------------
library(bbssr)
library(dplyr)
library(ggplot2)
library(tidyr)

## ----chisq_example------------------------------------------------------------
# Example: Chi-squared test
power_chisq <- BinaryPower(
  p1 = 0.6, p2 = 0.4, 
  N1 = 30, N2 = 30, 
  alpha = 0.025, 
  Test = 'Chisq'
)
print(paste("Chi-squared test power:", round(power_chisq, 3)))

## ----fisher_example-----------------------------------------------------------
# Example: Fisher exact test
power_fisher <- BinaryPower(
  p1 = 0.6, p2 = 0.4, 
  N1 = 30, N2 = 30, 
  alpha = 0.025, 
  Test = 'Fisher'
)
print(paste("Fisher exact test power:", round(power_fisher, 3)))

## ----fisher_midp_example------------------------------------------------------
# Example: Fisher mid-p test
power_midp <- BinaryPower(
  p1 = 0.6, p2 = 0.4, 
  N1 = 30, N2 = 30, 
  alpha = 0.025, 
  Test = 'Fisher-midP'
)
print(paste("Fisher mid-p test power:", round(power_midp, 3)))

## ----zpool_example------------------------------------------------------------
# Example: Z-pooled test
power_zpool <- BinaryPower(
  p1 = 0.6, p2 = 0.4, 
  N1 = 30, N2 = 30, 
  alpha = 0.025, 
  Test = 'Z-pool'
)
print(paste("Z-pooled test power:", round(power_zpool, 3)))

## ----boschloo_example---------------------------------------------------------
# Example: Boschloo test
power_boschloo <- BinaryPower(
  p1 = 0.6, p2 = 0.4, 
  N1 = 30, N2 = 30, 
  alpha = 0.025, 
  Test = 'Boschloo'
)
print(paste("Boschloo test power:", round(power_boschloo, 3)))

## ----test_comparison----------------------------------------------------------
# Compare all five tests
tests <- c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')
powers <- sapply(tests, function(test) {
  BinaryPower(p1 = 0.6, p2 = 0.4, N1 = 30, N2 = 30, alpha = 0.025, Test = test)
})

comparison_df <- data.frame(
  Test = tests,
  Power = round(powers, 4),
  Type = c("Asymptotic", "Conditional", "Conditional", "Unconditional", "Unconditional"),
  Conservatism = c("Moderate", "High", "Moderate", "Moderate", "Low")
)

print(comparison_df)

## ----bssr_detailed------------------------------------------------------------
# Detailed BSSR example with different approaches
bssr_results_list <- list()

# Restricted approach
bssr_results_list[["Restricted"]] <- BinaryPowerBSSR(
  asmd.p1 = 0.45, asmd.p2 = 0.09, 
  p = seq(0.1, 0.9, by = 0.1),
  Delta.A = 0.36, Delta.T = 0.36, 
  N1 = 24, N2 = 24, omega = 0.5, r = 1, 
  alpha = 0.025, tar.power = 0.8, 
  Test = 'Z-pool', 
  restricted = TRUE, weighted = FALSE
) %>% mutate(approach = "Restricted")

# Unrestricted approach
bssr_results_list[["Unrestricted"]] <- BinaryPowerBSSR(
  asmd.p1 = 0.45, asmd.p2 = 0.09, 
  p = seq(0.1, 0.9, by = 0.1),
  Delta.A = 0.36, Delta.T = 0.36, 
  N1 = 24, N2 = 24, omega = 0.5, r = 1, 
  alpha = 0.025, tar.power = 0.8, 
  Test = 'Z-pool', 
  restricted = FALSE, weighted = FALSE
) %>% mutate(approach = "Unrestricted")

# Weighted approach
bssr_results_list[["Weighted"]] <- BinaryPowerBSSR(
  asmd.p1 = 0.45, asmd.p2 = 0.09, 
  p = seq(0.1, 0.9, by = 0.1),
  Delta.A = 0.36, Delta.T = 0.36, 
  N1 = 24, N2 = 24, omega = 0.5, r = 1, 
  alpha = 0.025, tar.power = 0.8, 
  Test = 'Z-pool', 
  restricted = FALSE, weighted = TRUE
) %>% mutate(approach = "Weighted")

# Combine results
bssr_results <- do.call(rbind, bssr_results_list)

# Summary statistics
bssr_summary <- bssr_results %>%
  group_by(approach) %>%
  summarise(
    mean_power_bssr = mean(power.BSSR),
    mean_power_trad = mean(power.TRAD),
    min_power_bssr = min(power.BSSR),
    max_power_bssr = max(power.BSSR),
    .groups = 'drop'
  )

print(bssr_summary)

## ----power_comparison_plot, fig.width=8, fig.height=10, out.width="100%"------
# Create comprehensive power comparison with vertical layout
power_data <- bssr_results %>%
  select(approach, p, power.BSSR, power.TRAD) %>%
  pivot_longer(
    cols = c(power.BSSR, power.TRAD),
    names_to = "design_type",
    values_to = "power"
  ) %>%
  mutate(
    design_type = case_when(
      design_type == "power.BSSR" ~ "BSSR",
      design_type == "power.TRAD" ~ "Traditional"
    ),
    approach = factor(approach, levels = c("Restricted", "Unrestricted", "Weighted"))
  )

ggplot(power_data, aes(x = p, y = power, color = design_type)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~approach, ncol = 1, scales = "free_y") +  # Vertical layout
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray") +
  scale_color_manual(
    values = c("BSSR" = "#1F78B4", "Traditional" = "#E31A1C"),
    name = "Design Type"
  ) +
  scale_x_continuous(
    breaks = seq(0.2, 0.8, by = 0.2),
    labels = c("0.2", "0.4", "0.6", "0.8")
  ) +
  scale_y_continuous(
    breaks = seq(0.7, 1.0, by = 0.1),
    labels = c("0.7", "0.8", "0.9", "1.0")
  ) +
  labs(
    x = "Pooled Response Rate (θ)",
    y = "Power",
    title = "Power Comparison: Traditional vs BSSR Designs",
    subtitle = "Horizontal dashed line shows target power = 0.8"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 15)),
    strip.text = element_text(size = 12, face = "bold", margin = margin(t = 8, b = 8)),
    strip.background = element_rect(fill = "gray95", color = "gray80"),
    legend.position = "bottom",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.margin = margin(t = 10),
    axis.title.x = element_text(size = 11, margin = margin(t = 10)),
    axis.title.y = element_text(size = 11, margin = margin(r = 10)),
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray92", linewidth = 0.5),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

## ----sample_size_planning-----------------------------------------------------
# Sample size planning example
planning_scenarios <- expand.grid(
  p1 = c(0.4, 0.5, 0.6),
  p2 = c(0.2, 0.3),
  test = c('Fisher', 'Z-pool', 'Boschloo')
) %>%
  filter(p1 > p2)

# Calculate sample sizes for each scenario
sample_size_results <- list()
for(i in 1:nrow(planning_scenarios)) {
  result <- BinarySampleSize(
    p1 = planning_scenarios$p1[i], 
    p2 = planning_scenarios$p2[i], 
    r = 1, 
    alpha = 0.025, 
    tar.power = 0.8, 
    Test = planning_scenarios$test[i]
  )
  sample_size_results[[i]] <- result
}

# Combine results
final_results <- do.call(rbind, sample_size_results)
final_results <- final_results[, c("p1", "p2", "Test", "N1", "N2", "N", "Power")]

print(final_results)

## ----type1_error--------------------------------------------------------------
# Demonstrate Type I error control under null hypothesis
null_powers <- sapply(c('Fisher', 'Z-pool', 'Boschloo'), function(test) {
  BinaryPower(p1 = 0.3, p2 = 0.3, N1 = 30, N2 = 30, alpha = 0.025, Test = test)
})

names(null_powers) <- c('Fisher', 'Z-pool', 'Boschloo')
print("Type I error rates under null hypothesis:")
print(round(null_powers, 4))

## ----allocation_ratios--------------------------------------------------------
# Compare different allocation ratios
ratios <- c(1, 2, 3)
ratio_results <- sapply(ratios, function(r) {
  result <- BinarySampleSize(
    p1 = 0.5, p2 = 0.3, r = r, 
    alpha = 0.025, tar.power = 0.8, 
    Test = 'Boschloo'
  )
  c(N1 = result$N1, N2 = result$N2, N_total = result$N)
})

colnames(ratio_results) <- paste0("r=", ratios)
print("Sample sizes for different allocation ratios:")
print(ratio_results)

## ----sensitivity_analysis-----------------------------------------------------
# Sensitivity analysis for key parameters
sensitivity_data <- expand.grid(
  omega = c(0.3, 0.5, 0.7),
  alpha = c(0.01, 0.025, 0.05)
) %>%
  rowwise() %>%
  mutate(
    avg_power = mean(BinaryPowerBSSR(
      asmd.p1 = 0.45, asmd.p2 = 0.09, 
      p = seq(0.2, 0.8, by = 0.1),
      Delta.A = 0.36, Delta.T = 0.36, 
      N1 = 24, N2 = 24, omega = omega, r = 1, 
      alpha = alpha, tar.power = 0.8, 
      Test = 'Z-pool', 
      restricted = FALSE, weighted = FALSE
    )$power.BSSR)
  )

print("Sensitivity analysis results:")
print(sensitivity_data)

