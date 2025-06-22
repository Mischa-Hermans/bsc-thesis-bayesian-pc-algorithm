# ---- utils/prior_generation.R ----
# Prior Matrix Generation for PC Simulations
#
# Description:
#   Provides functions to generate prior knowledge matrices (fixedEdges and fixedGaps)
#   for both hard constraint and Bayesian approaches in PC algorithm simulations.
#   Supports tuning of true/false positive and negative rates for different prior qualities.
#
# Provides:
#   - sample_fixedGaps(): Samples a fixedGaps matrix with user-defined TN/FN rates
#   - sample_fixedEdges(): Samples a fixedEdges matrix with user-defined TP/FP rates
#   - generate_prior_list(): Returns a list of prior configurations for selected settings
#
# Dependencies:
#   - stats
#
# Author: Mischa Hermans

# Load required libraries
library(stats)

# Sample fixedGaps matrix given true/false negative rates
sample_fixedGaps <- function(A, pTN, pFN, symmetric = TRUE) {
  p <- nrow(A)
  fg <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      if (i == j) next
      if (A[i, j] == 0 && A[j, i] == 0 && runif(1) < pTN) fg[i, j] <- 1
      if ((A[i, j] != 0 || A[j, i] != 0) && runif(1) < pFN) fg[i, j] <- 1
    }
  }
  diag(fg) <- 0
  if (symmetric) fg <- 1 * ((fg + t(fg)) > 0)
  fg
}

# Sample fixedEdges matrix given true/false positive rates
sample_fixedEdges <- function(A, pTP, pFP, symmetric = TRUE) {
  fe <- matrix(0, nrow(A), ncol(A))
  for (i in 1:nrow(A)) {
    for (j in 1:ncol(A)) {
      if (i == j) next
      if (A[i, j] != 0 && runif(1) < pTP) fe[i, j] <- 1
      if (A[i, j] == 0 && runif(1) < pFP) fe[i, j] <- 1
    }
  }
  diag(fe) <- 0
  if (symmetric) fe <- 1 * ((fe + t(fe)) > 0)
  fe
}

# Generate prior configuration list for all supported settings
generate_prior_list <- function(A, pTP = 0.8, pFP = 0.2, settings_to_run) {
  fe_best <- sample_fixedEdges(A, 1, 0)
  fe_worst <- sample_fixedEdges(A, 0, 1)
  
  priors <- list(
    No_Background = list(),
    Hard_Best = list(fixedEdges = fe_best, fixedGaps = 1 * (fe_best == 0)),
    Hard_Realistic = list(
      fixedEdges = sample_fixedEdges(A, pTP / 2, pFP / 2),
      fixedGaps = sample_fixedGaps(A, pTP / 2, pFP / 2)
    ),
    Hard_Conservative = list(
      fixedEdges = sample_fixedEdges(A, 0.15, 0),
      fixedGaps = sample_fixedGaps(A, 0.15, 0)
    ),
    Hard_Worst = list(fixedEdges = fe_worst, fixedGaps = 1 * (fe_worst == 0)),
    
    Bayesian_Best = list(prior = sample_fixedEdges(A, 1, 0, FALSE)),
    Bayesian_Realistic = list(prior = sample_fixedEdges(A, pTP, pFP, FALSE)),
    Bayesian_Conservative = list(prior = sample_fixedEdges(A, 0.3, 0, FALSE)),
    Bayesian_Spike_Slab = list(prior = sample_fixedEdges(A, pTP, pFP, FALSE)),
    Bayesian_Worst = list(prior = sample_fixedEdges(A, 0, 1, FALSE))
  )
  
  priors[names(priors) %in% settings_to_run]
}
