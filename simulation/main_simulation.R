# ---- simulation/main_simulation.R ----
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
#   - models/stan_models.R
#   - utils/ci_tests.R
#   - utils/data_generation.R
#   - utils/evaluation_metrics.R
#   - utils/results_summary.R
#   - simulaton/simulation_functions.R
#   - pcalg, igraph, dplyr, rstan
#
# Author: Mischa Hermans

# Load required project files
source("models/stan_models.R")
source("utils/ci_tests.R")
source("utils/evaluation_metrics.R")
source("utils/results_summary.R")
source("simulation/simulation_functions.R")

# Load required libraries
library(pcalg)
library(igraph)
library(dplyr)
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
                     "Bayesian_Best", "Bayesian_Realistic", "Bayesian_Conservative", "Bayesian_Worst",
                     "Bayesian_Spike_Slab")

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

# Show results
generate_confusion_table(simulation_results$confusion)
generate_effects_table(simulation_results$effects)
