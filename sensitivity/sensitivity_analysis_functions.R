# ---- sensitivity_analysis_functions.R ----
# Sensitivity Analysis Functions for Causal Discovery Evaluation
#
# Description:
#   Provides core functions to run sensitivity analyses for the PC algorithm
#   under various prior knowledge settings and Bayesian CI tests. Each analysis
#   varies one parameter (e.g., pTP, pFP, mu, sigma, n, noise_sd) and evaluates
#   the impact on structure recovery and causal effect estimation.
#
# Functions:
#   - run_sensitivity_analysis():
#       Executes simulation runs for a single parameter grid, modifying only one
#       parameter at a time while keeping others fixed.
#       Returns confusion and effect metrics for each parameter value.
#
#   - aggregate_sensitivity_results():
#       Combines confusion and effect results over sensitivity analysis.
#
# Used by:
#   - run_sensitivity_analysis.R
#
# Depends on:
#   - simulate_methods() from simulation_functions.R
#   - igraph, pcalg, rstan
#
# Output:
#   - List of results per parameter value
#   - Aggregated summary tables (confusion and effects) per parameter
#
# Author: Mischa Hermans

# Runs a full sensitivity analysis by varying one parameter at a time.
run_sensitivity_analysis <- function(param_name, values,
                                     nruns = 10, alphas = c(0.2),
                                     p = 10, n = 150,
                                     discrete_idx = 1:3,
                                     A, true_graph,
                                     settings_to_run,
                                     var_types,
                                     stan_model,
                                     stan_model_spike_slab,
                                     base_params = list(pTP = 0.8, pFP = 0.2, mu = 1.5, sigma = 1, noise_sd = 1)) {
  results <- list()
  
  for (val in values) {
    cat(sprintf("\nRunning sensitivity for %s = %.2f\n", param_name, val))
    
    current_params <- base_params
    current_params[[param_name]] <- val
    current_n <- if (param_name == "n") val else n
    
    sim_result <- simulate_methods(
      nruns = nruns,
      alphas = alphas,
      p = p,
      n = current_n,
      A = A,
      true_graph = true_graph,
      discrete_idx = discrete_idx,
      pTP = current_params$pTP,
      pFP = current_params$pFP,
      mu = current_params$mu,
      sigma = current_params$sigma,
      noise_sd = current_params$noise_sd,
      settings_to_run = settings_to_run,
      var_types = var_types,
      stan_model = stan_model,
      stan_model_spike_slab = stan_model_spike_slab
    )
    
    results[[as.character(val)]] <- list(
      param_name = param_name,
      param_value = val,
      confusion = sim_result$confusion,
      effects = sim_result$effects
    )
  }
  
  results
}

# Aggregates confusion and effects results across all parameter values
aggregate_sensitivity_results <- function(sensitivity_results) {
  confusion_all <- do.call(rbind, lapply(sensitivity_results, function(res) {
    df <- res$confusion
    df$param_name <- res$param_name
    df$param_value <- res$param_value
    df
  }))
  
  effects_all <- do.call(rbind, lapply(sensitivity_results, function(res) {
    df <- res$effects
    df$param_name <- res$param_name
    df$param_value <- res$param_value
    df
  }))
  
  list(confusion = confusion_all, effects = effects_all)
}
