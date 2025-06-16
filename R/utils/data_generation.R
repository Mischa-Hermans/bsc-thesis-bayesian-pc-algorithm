# ---- R/data_generation.R ----
# Data Generation from DAG Structure
#
# Description:
#   Simulates a dataset from a known directed acyclic graph (DAG) with mixed data types.
#   Discrete variables are sampled as categorical, continuous variables follow linear Gaussian models.
#
# Provides:
#   - generate_data_from_dag(): Simulates n observations from DAG A, supporting both discrete and continuous variables.
#
# Dependencies:
#   - stats (rnorm, sample)
#
# Author: Mischa Hermans

generate_data_from_dag <- function(n, p, A, discrete_idx, noise_sd = 1) {
  X <- matrix(NA, nrow = n, ncol = p)
  
  for (j in 1:p) {
    parents <- which(A[, j] != 0)
    
    if (j %in% discrete_idx) {
      # Sample categorical variable with 3 levels
      X[, j] <- sample(1:3, size = n, replace = TRUE)
    } else {
      noise <- rnorm(n, mean = 0, sd = noise_sd)
      X[, j] <- if (length(parents) == 0) {
        rnorm(n)
      } else {
        X[, parents, drop = FALSE] %*% A[parents, j] + noise
      }
    }
  }
  
  colnames(X) <- colnames(A)
  return(X)
}
