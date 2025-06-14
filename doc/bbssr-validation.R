## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 6,
  fig.align = "center",
  warning = FALSE,
  message = FALSE
)

# Check if validation packages are available
has_exact <- requireNamespace("Exact", quietly = TRUE)
has_exact2x2 <- requireNamespace("exact2x2", quietly = TRUE)

if (!has_exact || !has_exact2x2) {
  knitr::opts_chunk$set(eval = FALSE)
  warning("Validation packages not available. Code will not be evaluated.")
}

## ----load_packages, eval=has_exact && has_exact2x2----------------------------
library(bbssr)
library(Exact)
library(exact2x2)
library(microbenchmark)
library(knitr)
library(ggplot2)
library(dplyr)

## ----load_bbssr_only, eval=!has_exact || !has_exact2x2------------------------
# library(bbssr)

## ----validation_parameters----------------------------------------------------
# Test parameters for validation
p1 <- c(0.5, 0.6, 0.7, 0.8)  # Response rates in treatment group
p2 <- c(0.2, 0.2, 0.2, 0.2)  # Response rates in control group
N1 <- 10                     # Sample size in treatment group
N2 <- 40                     # Sample size in control group
alpha <- 0.025              # One-sided significance level

## ----quick_verification-------------------------------------------------------
# Quick verification using bbssr only
cat("=== Quick Verification Results ===\n")
cat("Test Parameters:\n")
cat("p1 =", paste(p1, collapse = ", "), "\n")
cat("p2 =", paste(p2, collapse = ", "), "\n") 
cat("N1 =", N1, ", N2 =", N2, ", alpha =", alpha, "\n\n")

# Show all test results
tests <- c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')
for(test in tests) {
  powers <- BinaryPower(p1, p2, N1, N2, alpha, Test = test)
  cat(sprintf("%-12s: %s\n", test, paste(round(powers, 6), collapse = "  ")))
}

## ----chisq_validation, eval=has_exact && has_exact2x2-------------------------
# bbssr results
bbssr_chisq <- BinaryPower(p1, p2, N1, N2, alpha, Test = 'Chisq')
print("bbssr results:")
print(round(bbssr_chisq, 7))

# Exact package results for comparison
exact_chisq <- sapply(seq(p1), function(i) {
  power.exact.test(p1[i], p2[i], N1, N2, method = 'pearson chisq', 'greater', alpha)$power
})
print("Exact package results:")
print(round(exact_chisq, 7))

# Create comparison table
chisq_comparison <- data.frame(
  p1 = p1,
  p2 = p2,
  bbssr = round(bbssr_chisq, 7),
  Exact = round(exact_chisq, 7),
  Difference = round(abs(bbssr_chisq - exact_chisq), 10)
)

kable(chisq_comparison, caption = "Pearson Chi-squared Test: bbssr vs Exact Package")

## ----fisher_validation, eval=has_exact && has_exact2x2------------------------
# bbssr results
bbssr_fisher <- BinaryPower(p1, p2, N1, N2, alpha, Test = 'Fisher')

# Exact package results for comparison
exact_fisher <- sapply(seq(p1), function(i) {
  power.exact.test(p1[i], p2[i], N1, N2, method = 'fisher', 'greater', alpha)$power
})

# Create comparison table
fisher_comparison <- data.frame(
  p1 = p1,
  p2 = p2,
  bbssr = round(bbssr_fisher, 7),
  Exact = round(exact_fisher, 7),
  Difference = round(abs(bbssr_fisher - exact_fisher), 10)
)

kable(fisher_comparison, caption = "Fisher Exact Test: bbssr vs Exact Package")

## ----fisher_midp_validation, eval=has_exact && has_exact2x2-------------------
# bbssr results
bbssr_midp <- BinaryPower(p1, p2, N1, N2, alpha, Test = 'Fisher-midP')

# exact2x2 package results for comparison
exact2x2_midp <- sapply(seq(p1), function(i) {
  Power2x2(N1, N2, p1[i], p2[i], alpha, pvalFunc = 
             function(x1, n1, x2, n2) {
               fisher.exact(matrix(c(x1, n1 - x1, x2, n2 - x2), 2), 
                          alt = 'greater', midp = TRUE)$p.value
             }
  )
})

# Create comparison table
midp_comparison <- data.frame(
  p1 = p1,
  p2 = p2,
  bbssr = round(bbssr_midp, 7),
  exact2x2 = round(exact2x2_midp, 7),
  Difference = round(abs(bbssr_midp - exact2x2_midp), 10)
)

kable(midp_comparison, caption = "Fisher Mid-p Test: bbssr vs exact2x2 Package")

## ----zpool_validation, eval=has_exact && has_exact2x2-------------------------
# bbssr results
bbssr_zpool <- BinaryPower(p1, p2, N1, N2, alpha, Test = 'Z-pool')

# Exact package results for comparison
exact_zpool <- sapply(seq(p1), function(i) {
  power.exact.test(p1[i], p2[i], N1, N2, method = 'z-pooled', 'greater', alpha)$power
})

# Create comparison table
zpool_comparison <- data.frame(
  p1 = p1,
  p2 = p2,
  bbssr = round(bbssr_zpool, 7),
  Exact = round(exact_zpool, 7),
  Difference = round(abs(bbssr_zpool - exact_zpool), 10)
)

kable(zpool_comparison, caption = "Z-pooled Test: bbssr vs Exact Package")

## ----boschloo_validation, eval=has_exact && has_exact2x2----------------------
# bbssr results
bbssr_boschloo <- BinaryPower(p1, p2, N1, N2, alpha, Test = 'Boschloo')

# Exact package results for comparison
exact_boschloo <- sapply(seq(p1), function(i) {
  power.exact.test(p1[i], p2[i], N1, N2, method = 'boschloo', 'greater', alpha)$power
})

# Create comparison table
boschloo_comparison <- data.frame(
  p1 = p1,
  p2 = p2,
  bbssr = round(bbssr_boschloo, 7),
  Exact = round(exact_boschloo, 7),
  Difference = round(abs(bbssr_boschloo - exact_boschloo), 10)
)

kable(boschloo_comparison, caption = "Boschloo Test: bbssr vs Exact Package")

## ----sample_size_scenarios----------------------------------------------------
# Define test scenarios for sample size validation
scenarios <- list(
  list(p1 = 0.4, p2 = 0.2, r = 2, alpha = 0.025, tar.power = 0.8, Test = 'Chisq'),
  list(p1 = 0.5, p2 = 0.2, r = 3, alpha = 0.025, tar.power = 0.9, Test = 'Fisher'),
  list(p1 = 0.6, p2 = 0.3, r = 2, alpha = 0.025, tar.power = 0.9, Test = 'Fisher-midP'),
  list(p1 = 0.3, p2 = 0.2, r = 1, alpha = 0.025, tar.power = 0.8, Test = 'Z-pool'),
  list(p1 = 0.7, p2 = 0.2, r = 4, alpha = 0.025, tar.power = 0.8, Test = 'Boschloo')
)

## ----sample_size_validation, eval=has_exact && has_exact2x2-------------------
# Function to validate sample size calculation
validate_sample_size <- function(scenario) {
  # Calculate sample size using bbssr
  bbssr_result <- BinarySampleSize(
    p1 = scenario$p1, p2 = scenario$p2, r = scenario$r, 
    alpha = scenario$alpha, tar.power = scenario$tar.power, 
    Test = scenario$Test
  )
  
  # Verify power using external packages
  if (scenario$Test == 'Chisq') {
    external_power <- power.exact.test(
      scenario$p1, scenario$p2, bbssr_result$N1, bbssr_result$N2, 
      method = 'pearson chisq', 'greater', scenario$alpha
    )$power
  } else if (scenario$Test == 'Fisher') {
    external_power <- power.exact.test(
      scenario$p1, scenario$p2, bbssr_result$N1, bbssr_result$N2, 
      method = 'fisher', 'greater', scenario$alpha
    )$power
  } else if (scenario$Test == 'Fisher-midP') {
    external_power <- Power2x2(
      bbssr_result$N1, bbssr_result$N2, scenario$p1, scenario$p2, 
      alpha = scenario$alpha, 
      pvalFunc = function(x1, n1, x2, n2) {
        fisher.exact(matrix(c(x1, n1 - x1, x2, n2 - x2), 2), 
                    alt = 'greater', midp = TRUE)$p.value
      }
    )
  } else if (scenario$Test == 'Z-pool') {
    external_power <- power.exact.test(
      scenario$p1, scenario$p2, bbssr_result$N1, bbssr_result$N2, 
      method = 'z-pooled', 'greater', scenario$alpha
    )$power
  } else if (scenario$Test == 'Boschloo') {
    external_power <- power.exact.test(
      scenario$p1, scenario$p2, bbssr_result$N1, bbssr_result$N2, 
      method = 'boschloo', 'greater', scenario$alpha
    )$power
  }
  
  return(data.frame(
    Test = scenario$Test,
    p1 = scenario$p1,
    p2 = scenario$p2,
    r = scenario$r,
    Target_Power = scenario$tar.power,
    N1 = bbssr_result$N1,
    N2 = bbssr_result$N2,
    N_Total = bbssr_result$N,
    bbssr_Power = round(bbssr_result$Power, 6),
    External_Power = round(external_power, 6),
    Power_Difference = round(abs(bbssr_result$Power - external_power), 8)
  ))
}

# Run validation for all scenarios
validation_results <- do.call(rbind, lapply(scenarios, validate_sample_size))
kable(validation_results, caption = "Sample Size Validation: bbssr vs External Packages")

## ----performance_setup--------------------------------------------------------
# Parameters for performance comparison (smaller for exact2x2)
perf_p1 <- 0.6
perf_p2 <- 0.3
perf_N1 <- 8   # Smaller sample size for exact2x2 comparison
perf_N2 <- 12
perf_alpha <- 0.025

## ----power_performance, eval=has_exact && has_exact2x2------------------------
# Benchmark BinaryPower function
power_benchmark <- microbenchmark(
  bbssr_chisq = BinaryPower(perf_p1, perf_p2, perf_N1, perf_N2, perf_alpha, Test = 'Chisq'),
  exact_chisq = power.exact.test(perf_p1, perf_p2, perf_N1, perf_N2, 
                                method = 'pearson chisq', 'greater', perf_alpha),
  
  bbssr_fisher = BinaryPower(perf_p1, perf_p2, perf_N1, perf_N2, perf_alpha, Test = 'Fisher'),
  exact_fisher = power.exact.test(perf_p1, perf_p2, perf_N1, perf_N2, 
                                 method = 'fisher', 'greater', perf_alpha),
  
  bbssr_midp = BinaryPower(perf_p1, perf_p2, perf_N1, perf_N2, perf_alpha, Test = 'Fisher-midP'),
  exact2x2_midp = Power2x2(perf_N1, perf_N2, perf_p1, perf_p2, alpha = perf_alpha, 
                          pvalFunc = function(x1, n1, x2, n2) {
                            fisher.exact(matrix(c(x1, n1 - x1, x2, n2 - x2), 2), 
                                       alt = 'greater', midp = TRUE)$p.value
                          }),
  
  bbssr_zpool = BinaryPower(perf_p1, perf_p2, perf_N1, perf_N2, perf_alpha, Test = 'Z-pool'),
  exact_zpool = power.exact.test(perf_p1, perf_p2, perf_N1, perf_N2, 
                                method = 'z-pooled', 'greater', perf_alpha),
  
  bbssr_boschloo = BinaryPower(perf_p1, perf_p2, perf_N1, perf_N2, perf_alpha, Test = 'Boschloo'),
  exact_boschloo = power.exact.test(perf_p1, perf_p2, perf_N1, perf_N2, 
                                   method = 'boschloo', 'greater', perf_alpha),
  
  times = 50
)

# Display benchmark results
print(power_benchmark)

## ----performance_plot, eval=has_exact && has_exact2x2, fig.width=10, fig.height=6----
# Create performance comparison plot
benchmark_summary <- summary(power_benchmark)
benchmark_df <- data.frame(
  Test = c("Chi-squared", "Chi-squared", "Fisher", "Fisher", 
           "Fisher-midP", "Fisher-midP", "Z-pooled", "Z-pooled", 
           "Boschloo", "Boschloo"),
  Package = rep(c("bbssr", "External"), 5),
  Median_Time = benchmark_summary$median / 1000  # Convert to microseconds
)

ggplot(benchmark_df, aes(x = Test, y = Median_Time, fill = Package)) +
  geom_col(position = "dodge") +
  scale_y_log10() +
  labs(
    title = "Performance Comparison: bbssr vs External Packages",
    subtitle = "Median execution time (microseconds, log scale)",
    x = "Statistical Test",
    y = "Median Time (microseconds)",
    fill = "Package"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

## ----speed_improvement, eval=has_exact && has_exact2x2------------------------
# Calculate speed improvement ratios
speed_ratios <- data.frame(
  Test = c("Chi-squared", "Fisher", "Fisher-midP", "Z-pooled", "Boschloo"),
  bbssr_time = benchmark_summary$median[c(1, 3, 5, 7, 9)],
  external_time = benchmark_summary$median[c(2, 4, 6, 8, 10)],
  stringsAsFactors = FALSE
) %>%
  mutate(
    Speed_Improvement = round(external_time / bbssr_time, 2),
    bbssr_ms = round(bbssr_time / 1000, 3),
    external_ms = round(external_time / 1000, 3)
  ) %>%
  select(Test, bbssr_ms, external_ms, Speed_Improvement)

kable(speed_ratios, 
      col.names = c("Test", "bbssr (ms)", "External (ms)", "Speed Improvement (x)"),
      caption = "Speed Improvement: bbssr vs External Packages")

## ----validation_summary, eval=has_exact && has_exact2x2-----------------------
avg_improvement <- mean(speed_ratios$Speed_Improvement)
max_improvement <- max(speed_ratios$Speed_Improvement)
min_improvement <- min(speed_ratios$Speed_Improvement)

cat("Performance Summary:\n")
cat(sprintf("- Average speed improvement: %.1fx faster\n", avg_improvement))
cat(sprintf("- Maximum speed improvement: %.1fx faster\n", max_improvement))
cat(sprintf("- Minimum speed improvement: %.1fx faster\n", min_improvement))

## ----session_info-------------------------------------------------------------
sessionInfo()

