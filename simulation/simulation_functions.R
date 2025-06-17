# ---- simulation_functions.R ----
# Simulation Runner for Causal Discovery under Prior Knowledge
#
# Description:
#   Contains the core simulation loop for running the PC algorithm across different
#   prior knowledge configurations and alpha levels. Handles both hard constraints
#   and Bayesian CI testing using compiled Stan models.
#
# Provides:
#   - run_simulation(): Runs a single simulation across multiple settings and alphas
#   - simulate_methods(): Aggregates results across multiple simulation runs
#
# Dependencies:
#   - igraph
#   - pcalg
#   - stats (cor)
#   - data_generation.R
#   - prior_generation.R
#   - evaluation_metrics.R
#   - ci_tests.R (for bayesCItest_symmetric and bayesCItest_spike_slab_symmetric)
#
# Author: Mischa Hermans

source("prior_generation.R")
source("data_generation.R")

# Executes a single simulation run for multiple alpha levels and prior settings
run_simulation <- function(run, alphas, p, n, A, true_graph, discrete_idx, var_types, 
                           pTP = 0.8, pFP = 0.2, mu = 1.5, sigma = 1, noise_sd = 1,
                           settings_to_run, stan_model, stan_model_spike_slab) {
  all_confusion <- data.frame()
  all_effects <- data.frame()
  
  # Simulate dataset X according to true dag A
  X <- generate_data_from_dag(n, p, A, discrete_idx, noise_sd)
  
  # Generate the priors
  priors <- generate_prior_list(A, pTP, pFP, settings_to_run)

  # Loop over significance levels and settings
  for (alpha in alphas) {
    for (setting in names(priors)) {
      
      if (setting == "No_Background") {
        suffStat <- list(C = cor(X), n = n)
        pc_fit <- pc(suffStat, indepTest = gaussCItest, alpha = alpha, labels = colnames(A))
      } 
      
      if (grepl("^Hard", setting)) {
        suffStat <- list(C = cor(X), n = n)
        pc_fit <- pc(suffStat, indepTest = gaussCItest, alpha = alpha, labels = colnames(A),
                     fixedEdges = priors[[setting]]$fixedEdges, fixedGaps = priors[[setting]]$fixedGaps)
      }
      
      if (grepl("^Bayesian", setting)) {
        prior_mask <- priors[[setting]]$prior
        mu_mat <- matrix(0, p, p)
        sigma_mat <- matrix(5, p, p)
        
        for (i in 1:p) for (j in 1:p) {
          if (prior_mask[i, j] == 1) {
            mu_mat[i,j] <- mu
            sigma_mat[i,j] <- sigma
          }
        }
        
        cache <- new.env(hash = TRUE, parent = emptyenv())
        
        if (setting == "Bayesian_Spike_Slab") {
          suffStat <- list(data = X,
                           prior_mu = mu_mat,
                           prior_sigma = sigma_mat,
                           prior_mask = prior_mask,
                           stan_model_spike_slab = stan_model_spike_slab,
                           alpha = alpha,
                           cache = cache)
          pc_fit <- pc(suffStat, indepTest = bayesCItest_spike_slab_symmetric, alpha = alpha, labels = colnames(A))
        } else {
          suffStat <- list(data = X,
                           prior_mu = mu_mat,
                           prior_sigma = sigma_mat,
                           stan_model = stan_model,
                           alpha = alpha,
                           cache = cache)

          pc_fit <- pc(suffStat, indepTest = bayesCItest_symmetric, alpha = alpha, labels = colnames(A))
        }
      }
      
      # Store results
      g_est <- graph_from_graphnel(pc_fit@graph)
      adj_est <- as.matrix(as_adjacency_matrix(g_est, sparse = FALSE))
      
      all_confusion <- rbind(all_confusion,
                             cbind(run = run, alpha = alpha, setting = setting,
                                   track_confusion(true_graph, g_est, var_types)))
      
      all_effects <- rbind(all_effects,
                           cbind(run = run, alpha = alpha, setting = setting,
                                 track_effects(X, adj_est, A, var_types)))
      
      confusion_part <- cbind(run = run, alpha = alpha, setting = setting,
                              track_confusion(true_graph, g_est, var_types))
      
      all_confusion <- rbind(all_confusion, confusion_part)
      
      confusion_all <- subset(confusion_part, category == "all")
      
      cat(sprintf(
        "  Setting: %s (alpha = %.2f) => TP = %d, FP = %d, FN = %d, TN = %d\n",
        setting, alpha,
        confusion_all$TP,
        confusion_all$FP,
        confusion_all$FN,
        confusion_all$TN
      ))
    }
  }
  
  list(confusion = all_confusion, effects = all_effects)
}

# Function to run multiple simulations
simulate_methods <- function(nruns, alphas, p, n, A, true_graph, discrete_idx,
                             pTP = 0.8, pFP = 0.2, mu = 1.5, sigma = 1, noise_sd = 1,
                             settings_to_run = NULL, var_types,
                             stan_model, stan_model_spike_slab) {
  
  all_confusion <- data.frame()
  all_effects <- data.frame()
  
  for (run in seq_len(nruns)) {
    cat(sprintf("\n--- RUN %d ---\n", run))
    result <- run_simulation(run, alphas, p, n, A, true_graph, discrete_idx, var_types,
                             pTP, pFP, mu, sigma, noise_sd, settings_to_run,
                             stan_model = stan_model, stan_model_spike_slab = stan_model_spike_slab)
    all_confusion <- rbind(all_confusion, result$confusion)
    all_effects   <- rbind(all_effects, result$effects)
  }
  
  list(confusion = all_confusion, effects = all_effects)
}
