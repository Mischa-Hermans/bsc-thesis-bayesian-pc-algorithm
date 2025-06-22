# ---- utils/summarize_results.R ----
# Summary Tables for Structure Recovery and Effect Estimation
#
# Description:
#   Generates formatted summary tables for confusion matrices and
#   causal effect estimation (MSE mean and MSE variance) by category.
#   Intended to be called after simulations have been run.
#
# Provides:
#   - summarize_confusion_table(): Summary of TP/FP/FN/TN by group
#   - summarize_effect_table(): Summary of MSE mean and MSE variance
#
# Dependencies:
#   - dplyr
#   - tidyr
#   - kableExtra
#
# Author: Mischa Hermans

# Load required libraries
library(dplyr)
library(tidyr)
library(kableExtra)

generate_confusion_table <- function(confusion_df) {
  confusion_df %>%
    group_by(alpha, setting, category) %>%
    summarize(across(c(TP, FP, FN, TN), mean), .groups = "drop") %>%
    complete(category = c("all", "discrete", "cts_indep", "cts_dep", "mixed")) %>%
    kable(digits = 1, caption = "Confusion Matrices by Alpha, Method, and Category") %>%
    kable_styling(full_width = FALSE)
}

generate_effects_table <- function(effects_df) {
  effects_df %>%
    group_by(alpha, setting, category) %>%
    summarize(
      mse_mean = mean(mse, na.rm = TRUE),
      mse_variance = var(mse, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    complete(category = c("all", "discrete", "cts_indep", "cts_dep", "mixed")) %>%
    kable(digits = 6, caption = "Effect Estimation: MSE Mean and Variance") %>%
    kable_styling(full_width = FALSE)
}

