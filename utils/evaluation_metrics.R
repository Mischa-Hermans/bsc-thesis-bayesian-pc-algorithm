# ---- utils/evaluation_metrics.R ----
# Evaluation Metrics for Structure and Effect Estimation
#
# Description:
#   Provides functions to compute confusion matrices and evaluate
#   estimation accuracy of causal effects across variable types.
#
# Provides:
#   - compute_confusion_matrix(): Basic TP/FP/FN/TN statistics
#   - compute_confusion_matrix_from_edges(): Confusion by edge type
#   - track_confusion(): Type-specific structure recovery evaluation
#   - track_effects(): MSE/variance for estimated regression effects
#
# Dependencies:
#   - igraph
#   - stats (lm)
#
# Author: Mischa Hermans

# Load required libraries
library(igraph)
library(stats)

# Compute basic confusion matrix (TP, FP, FN, TN) comparing true and estimated graphs
compute_confusion_matrix <- function(true_graph, est_graph) {
  true_adj <- as.matrix(as_adjacency_matrix(true_graph, sparse = FALSE))
  est_adj  <- as.matrix(as_adjacency_matrix(est_graph, sparse = FALSE))
  diag(true_adj) <- 0
  diag(est_adj)  <- 0
  
  # Flatten adjacency matrices into vectors
  true_vec <- as.vector(true_adj)
  est_vec  <- as.vector(est_adj)
  
  # Count TP/FP/FN/TN
  TP <- sum(true_vec == 1 & est_vec == 1)
  FP <- sum(true_vec == 0 & est_vec == 1)
  FN <- sum(true_vec == 1 & est_vec == 0)
  TN <- sum(true_vec == 0 & est_vec == 0)
  
  list(TP = TP, FP = FP, FN = FN, TN = TN)
}

# Compute confusion matrix using a mask to subset entries of adjacency matrices
compute_confusion_matrix_from_edges <- function(true_mat, est_mat, edge_mask) {
  true_vec <- as.vector(true_mat)[as.vector(edge_mask) == 1]
  est_vec  <- as.vector(est_mat)[as.vector(edge_mask) == 1]
  
  TP <- sum(true_vec == 1 & est_vec == 1)
  FP <- sum(true_vec == 0 & est_vec == 1)
  FN <- sum(true_vec == 1 & est_vec == 0)
  TN <- sum(true_vec == 0 & est_vec == 0)
  
  list(TP = TP, FP = FP, FN = FN, TN = TN)
}

# Compute confusion matrix statistics separately for variable type combinations
track_confusion <- function(true_graph, est_graph, types) {
  categories <- c("all", "discrete", "cts_indep", "cts_dep", "mixed")
  res <- data.frame()
  
  true_mat <- as.matrix(as_adjacency_matrix(true_graph, sparse = FALSE))
  est_mat <- as.matrix(as_adjacency_matrix(est_graph, sparse = FALSE))
  p <- ncol(true_mat)
  
  # Build binary masks for edge type combinations
  mask_discrete <- outer(types, types, function(a, b) a == "discrete" & b == "discrete")
  mask_cts_indep <- outer(types, types, function(a, b) a == "cts_indep" & b == "cts_indep")
  mask_cts_dep <- outer(types, types, function(a, b) a == "cts_dep" & b == "cts_dep")
  mask_mixed <- outer(types, types, function(a, b) {
    (a == "discrete" & (b == "cts_indep" | b == "cts_dep")) |
      (b == "discrete" & (a == "cts_indep" | a == "cts_dep")) |
      (a == "cts_indep" & b == "cts_dep") |
      (a == "cts_dep" & b == "cts_indep")
  })
  
  masks <- list(
    discrete = mask_discrete,
    cts_indep = mask_cts_indep,
    cts_dep = mask_cts_dep,
    mixed = mask_mixed
  )
  
  # Compute confusion matrix for each category
  for (cat in categories) {
    if (cat == "all") {
      cm <- compute_confusion_matrix(true_graph, est_graph)
    } else {
      cm <- compute_confusion_matrix_from_edges(true_mat, est_mat, masks[[cat]])
    }
    res <- rbind(res, data.frame(category = cat, cm))
  }
  res
}

# Evaluate MSE and variance of estimated regression effects for each variable type
track_effects <- function(X, adj_est, A, types) {
  categories <- c("all", unique(types), "mixed")
  res <- data.frame()
  
  for (cat in categories) {
    if (cat == "all") {
      idx <- 1:ncol(X)
    } else if (cat == "mixed") {
      # Identify variables with mixed-type parent-child relationships
      idx_pairs <- which(outer(types, types, FUN = "!=") & !diag(ncol(X)), arr.ind = TRUE)
      idx_list <- unique(c(idx_pairs[, 1], idx_pairs[, 2]))
      idx <- sort(unique(idx_list))
    } else {
      idx <- which(types == cat)
    }
    
    mse_list <- c()
    variance_list <- c()
    for (j in idx) {
      parents <- which(adj_est[, j] == 1)
      parents <- intersect(parents, idx)  # restrict to same category
      if (length(parents) > 0) {
        fit <- lm(X[, j] ~ X[, parents, drop = FALSE])
        est_effects <- coef(fit)[-1]
        true_effects <- A[parents, j]
        
        if (length(est_effects) > 0 && length(true_effects) == length(est_effects)) {
          mse_list <- c(mse_list, mean((est_effects - true_effects)^2, na.rm = TRUE))
          variance_list <- c(variance_list, var(est_effects, na.rm = TRUE))
        }
      }
    }
    
    res <- rbind(res, data.frame(
      category = cat,
      mse = if (length(mse_list) > 0) mean(mse_list, na.rm = TRUE) else NA,
      variance = if (length(variance_list) > 0) mean(variance_list, na.rm = TRUE) else NA
    ))
  }
  return(res)
}
