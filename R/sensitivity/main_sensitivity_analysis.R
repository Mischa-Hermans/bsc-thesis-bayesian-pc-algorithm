# ---- run_sensitivity_analysis.R ----
# Sensitivity Analysis Main Script for Causal Discovery Methods
#
# Description:
#   This script runs a full sensitivity analysis by varying key simulation parameters
#   (e.g., prior strength, noise levels, sample size) and evaluating their effect on
#   structure recovery and causal effect estimation.
#
#   Results across all parameters and settings are aggregated into summary tables.
#   The script outputs a combined structure recovery table (TP/FP/FN/TN) and
#   a causal effect estimation table (MSE, variance), filtered to category == "all".
#
# Features:
#   - Iterates over a grid of sensitivity parameters (e.g., pTP, pFP, mu, sigma, n)
#   - Evaluates multiple configurations of prior knowledge (e.g., No_Background, Hard, Bayesian)
#   - Compiles and uses Stan-based Bayesian CI tests
#   - Outputs results tables
#
# Depends on:
#   - R/stan_models.R             (defines and compiles Stan models)
#   - R/ci_tests.R                (Bayesian conditional independence tests)
#   - R/evaluation_metrics.R     (TP/FP/FN/TN, MSE, etc.)
#   - R/sensitivity_analysis_functions.R (core simulation loop and aggregation)
#
# Author: Mischa Hermans

# Load modules
source("R/stan_models.R")
source("R/ci_tests.R")
source("R/evaluation_metrics.R")
source("R/sensitivity_analysis_functions.R")

# Load required packages
library(pcalg)
library(igraph)
library(ggplot2)
library(purrr)
library(rstan)
library(dplyr)
library(tidyr)
library(knitr)

# Compile Stan models and define conditional independence tests
stan_models <- get_stan_models()
stan_model <- stan_models$standard
stan_model_spike_slab <- stan_models$spike_slab

ci_tests <- define_symmetric_tests(stan_model, stan_model_spike_slab)
bayesCItest_symmetric <- ci_tests$bayesCItest_symmetric
bayesCItest_spike_slab_symmetric <- ci_tests$bayesCItest_spike_slab_symmetric

# Define parameter grid for sensitivity analysis
sensitivity_specs <- list(
  pTP      = seq(0, 1, 0.1),
  pFP      = seq(0, 1, 0.1),
  sigma    = seq(0.2, 3, 0.3),
  n        = seq(100, 1000, 100),
  noise_sd = c(0.2, 0.5, 1, 1.5),
  mu       = c(0.5, 1, 1.5, 2)
)

# Define settings to run
settings_to_run <- c("No_Background", "Hard_Realistic", "Bayesian_Realistic", "Bayesian_Spike_Slab")

# Variable types
var_types <- rep("mixed", p)
var_types[discrete_idx] <- "discrete"
var_types[setdiff(1:p, c(discrete_idx, 10))] <- "cts_indep"
var_types[10] <- "cts_dep"

# Aggregate storage
confusion_all_combined <- data.frame()
effects_all_combined <- data.frame()

# Run over each parameter
for (param in names(sensitivity_specs)) {
  cat(sprintf("\n--- Running sensitivity for parameter: %s ---\n", param))
  
  results <- run_sensitivity_analysis(
    param_name = param,
    values = sensitivity_specs[[param]],
    nruns = 1,
    alphas = c(0.2),
    p = p,
    n = n,
    discrete_idx = discrete_idx,
    A = A,
    true_graph = true_graph,
    settings_to_run = settings_to_run,
    var_types = var_types,
    stan_model = stan_model,
    stan_model_spike_slab = stan_model_spike_slab
  )
  
  aggregated <- aggregate_sensitivity_results(results)
  
  confusion_all_combined <- rbind(
    confusion_all_combined,
    aggregated$confusion %>%
      filter(category == "all") %>%
      group_by(param_name, param_value, setting) %>%
      summarize(across(c(TP, FP, FN, TN), mean), .groups = "drop")
  )
  
  effects_all_combined <- rbind(
    effects_all_combined,
    aggregated$effects %>%
      filter(category == "all") %>%
      group_by(param_name, param_value, setting) %>%
      summarize(
        mse = mean(mse, na.rm = TRUE),
        variance = mean(variance, na.rm = TRUE),
        .groups = "drop"
      )
  )
}

# Print full structure recovery table
cat("\n=== Combined Structure Recovery Summary (category = 'all') ===\n")
print(kable(confusion_all_combined, digits = 2, caption = "Structure Recovery (All)"))

# Print full effect estimation table
cat("\n=== Combined Causal Effect Estimation Summary (category = 'all') ===\n")
print(kable(effects_all_combined, digits = 3, caption = "Causal Effect Estimation (All)"))
