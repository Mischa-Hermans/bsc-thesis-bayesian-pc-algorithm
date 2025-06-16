# ---- run_simulation.R ----
# Main Script to Run Causal Discovery Simulations and Summarize Results
#
# Description:
#   Executes the full simulation process for evaluating PC algorithm variants 
#   under different background knowledge settings (hard constraints and Bayesian CI tests).
#   Summarizes structure recovery and causal effect estimation across settings.
#
# Provides:
#   - Compilation of Stan models
#   - Definition of symmetric CI tests
#   - Specification of DAG, variable types, and simulation settings
#   - Simulation runs across multiple alpha levels and settings
#   - Aggregated performance summaries (TP/FP/FN/TN, MSE, Variance)
#
# Dependencies:
#   - R/stan_models.R
#   - R/ci_tests.R
#   - R/data_generation.R
#   - R/evaluation_metrics.R
#   - R/simulation_functions.R
#   - igraph
#   - pcalg
#   - rstan
#   - dplyr, tidyr, purrr, ggplot2, openxlsx, kableExtra
#
# Author: Mischa Hermans

# Load required project files
source("R/stan_models.R")
source("R/ci_tests.R")
source("R/evaluation_metrics.R")
source("R/simulation_functions.R")

# Load required libraries
library(pcalg)
library(igraph)
library(openxlsx)
library(ggplot2)
library(purrr)
library(tidyr)
library(dplyr)
library(kableExtra)
library(rstan)

# Compile Stan models and define conditional independence tests
stan_models <- get_stan_models()
stan_model <- stan_models$standard
stan_model_spike_slab <- stan_models$spike_slab

ci_tests <- define_symmetric_tests(stan_model, stan_model_spike_slab)
bayesCItest_symmetric <- ci_tests$bayesCItest_symmetric
bayesCItest_spike_slab_symmetric <- ci_tests$bayesCItest_spike_slab_symmetric

# Simulation parameters
nruns <- 100                 # Number of simulation repetitions
alphas <- c(0.2, 0.5, 0.8)   # Alpha levels for PC algorithm
p <- 10                      # Number of variables
n <- 150                     # Sample size
discrete_idx <- 1:3          # Indices of discrete variables

# Create the true DAG adjacency matrix
A <- matrix(0, p, p)
for (i in 1:9) A[i, 10] <- 1.5     # X1 to X9 cause Y
A[7, 8] <- 1.5                     # X7 → X8
A[7, 9] <- 1.5                     # X7 → X9
colnames(A) <- rownames(A) <- c(paste0("X", 1:9), "Y")

# Construct the true graph for evaluation (igraph object)
true_graph <- {
  g <- matrix(0, p, p)
  for (i in 1:9) g[i, 10] <- 1
  g[7, 8] <- 1
  g[7, 9] <- 1
  graph_from_adjacency_matrix(g, mode = "directed")
}

# Specify the variable types
var_types <- rep("mixed", p)
var_types[discrete_idx] <- "discrete"
var_types[setdiff(1:p, c(discrete_idx, 10))] <- "cts_indep"
var_types[10] <- "cts_dep"

# Specify which settings to run
settings_to_run <- c("No_Background", "Hard_Best", "Hard_Realistic", "Hard_Conservative", "Hard_Worst",
                     "Bayesian_Best", "Bayesian_Realistic", "Bayesian_Conservative", "Bayesian_Worst")

# Run the simulation
simulation_results <- simulate_methods(
  nruns = nruns,
  alphas = alphas,
  p = p,
  n = n,
  A = A,
  true_graph = true_graph,
  discrete_idx = discrete_idx,
  pTP = 0.8,
  pFP = 0.2,
  mu = 1.5,
  sigma = 1,
  noise_sd = 1,
  settings_to_run = settings_to_run,
  var_types = var_types,
  stan_model = stan_model,
  stan_model_spike_slab = stan_model_spike_slab
)

# Summarize structure recovery performance
confusion_summary <- simulation_results$confusion %>%
  group_by(alpha, setting, category) %>%
  summarize(across(c(TP, FP, FN, TN), mean), .groups = "drop") %>%
  complete(category = c("all", "discrete", "cts_indep", "cts_dep", "mixed"))

kable(confusion_summary, digits = 1, caption = "Confusion Matrices by Alpha, Method, and Category") %>%
  kable_styling(full_width = FALSE)

# Summarize causal effect estimation performance
effects_summary <- simulation_results$effects %>%
  group_by(alpha, setting, category) %>%
  summarize(
    mse = mean(mse, na.rm = TRUE),
    variance = mean(variance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(category = c("all", "discrete", "cts_indep", "cts_dep", "mixed"))

kable(effects_summary, digits = 3, caption = "MSE and Variance by Alpha, Method, and Category") %>%
  kable_styling(full_width = FALSE)
