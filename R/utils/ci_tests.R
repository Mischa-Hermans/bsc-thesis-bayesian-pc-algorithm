# ---- R/ci_tests.R ----
# Conditional Independence Tests for PC Algorithm
#
# Description:
#   Defines conditional independence tests based on Bayesian linear regression,
#   including standard Gaussian priors and spike-and-slab priors with inclusion weights.
#   Wraps the tests in symmetric wrappers and supports caching for efficiency.
#
# Provides:
#   - bayesCItest(): Bayesian CI test with standard normal priors
#   - bayesCItest_spike_slab(): Bayesian CI test with spike-and-slab priors
#   - make_symmetric_test(): Symmetrizes any CI test for use in PC algorithm
#   - define_symmetric_tests(): Returns all symmetric CI test variants with Stan model injection
#
# Dependencies:
#   - rstan
#   - stats
#
# Author: Mischa Hermans

# Performs a Bayesian conditional independence test using standard normal priors.
# Returns a p-value for the coefficient of x in the regression of y on x and S.
bayesCItest <- function(x, y, S, suffStat, epsilon = 0.1) {
  
  data <- suffStat$data
  prior_mu <- suffStat$prior_mu
  prior_sigma <- suffStat$prior_sigma
  n <- nrow(data)
  yvec <- data[, y]
  
  # Build regression design matrix using x and S as predictors
  predictors <- if (length(S)) c(x, S) else x
  Xmat <- data[, predictors, drop = FALSE]
  x_index <- which(predictors == x)  # Identify the index of x in the predictor list
  stopifnot(length(x_index) == 1)
  
  # Extract prior means and standard deviations for predictors
  mu_vec <- sapply(predictors, function(pred) prior_mu[pred, y])
  sigma_vec <- sapply(predictors, function(pred) prior_sigma[pred, y])
  mu_used <- mu_vec[x_index]
  
  # Format data for Stan
  stan_data <- list(
    N = n,
    K = ncol(Xmat),
    X = Xmat,
    y = yvec,
    mu_beta = as.array(mu_vec),
    sigma_beta = as.array(sigma_vec)
  )
 
  # Run the Stan model with prior information
  fit <- tryCatch({
    suppressMessages(suppressWarnings(
      sampling(suffStat$stan_model, data = stan_data, iter = 500, warmup = 250,
               chains = 1, refresh = 0, seed = 123, control = list(adapt_delta = 0.9))
    ))
  }, error = function(e) return(NULL))
  
  # If fitting fails, return NA p-value
  if (is.null(fit)) return(list(p = NA, mu = mu_used))
  
  # Extract posterior samples for the coefficient of x
  beta_post <- rstan::extract(fit)$beta[, x_index]
  if (length(beta_post) == 0 || all(is.na(beta_post))) return(list(p = NA, mu = mu_used))
  
  # Compute the two-sided bayesian p-value
  p_val <- 2 * min(mean(beta_post <= 0), mean(beta_post > 0))
  
  return(list(p = p_val, mu = mu_used))
}

# Performs a spike-and-slab Bayesian CI test with soft inclusion priors.
# Prior inclusion probabilities are derived from prior_mask.
bayesCItest_spike_slab <- function(x, y, S, suffStat, epsilon = 0.1) {
  data <- suffStat$data
  yvec <- data[, y]
  
  # Build regression design matrix using x and S
  predictors <- if (length(S)) c(x, S) else x
  Xmat <- data[, predictors, drop = FALSE]
  x_index <- which(predictors == x)
  
  # Build prior parameters for each predictor
  mu_vec <- sapply(predictors, function(pred) suffStat$prior_mu[pred, y])
  sigma_vec <- sapply(predictors, function(pred) suffStat$prior_sigma[pred, y])
  incl_vec <- sapply(predictors, function(pred) {
    if (suffStat$prior_mask[pred, y] == 1) 0.9 else 0.1
  })
  
  # Format input for the spike-and-slab Stan model
  stan_data <- list(
    N = nrow(data),
    K = length(predictors),
    X = Xmat,
    y = yvec,
    mu_beta = as.array(mu_vec),
    sigma_beta = as.array(sigma_vec),
    prior_inclusion = as.array(incl_vec)
  )
  
  # Run Stan model with inclusion weights
  fit <- tryCatch({
    suppressMessages(suppressWarnings(
      sampling(suffStat$stan_model_spike_slab, data = stan_data, iter = 500, warmup = 250,
               chains = 1, refresh = 0, seed = 123, control = list(adapt_delta = 0.95))
    ))
  }, error = function(e) return(NULL))
  
  # Handle failed fits gracefully
  if (is.null(fit)) return(list(p = NA, mu = NA))
  
  beta_post <- rstan::extract(fit)$beta[, x_index]
  if (length(beta_post) == 0 || all(is.na(beta_post))) return(list(p = NA, mu = NA))
  
  # Return two-sided p-value
  list(p = 2 * min(mean(beta_post <= 0), mean(beta_post > 0)), mu = NA)
}

# Wraps any CI test in a symmetric wrapper for use in the PC algorithm.
make_symmetric_test <- function(ci_func) {
  function(x, y, S, suffStat, alpha = 0.1) {
    cache <- suffStat$cache
    
    # Create a unique key based on variable pair and conditioning set
    S_key <- if (length(S) == 0) "∅" else paste(sort(S), collapse = ",")
    pair_key <- paste(sort(c(x, y)), collapse = "-")
    key <- paste0(pair_key, "|", S_key)
    
    # Check if result already cached
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache))
    }
    
    # Run the test in both directions and check if both p-values > alpha
    result1 <- ci_func(x, y, S, suffStat)
    result2 <- ci_func(y, x, S, suffStat)
    
    p1 <- result1$p
    p2 <- result2$p
    decision <- if (!is.na(p1) && !is.na(p2) && p1 > alpha && p2 > alpha) 1 else 0
    
    assign(key, decision, envir = cache)
    return(decision)
  }
}

# Returns a list of symmetric CI test functions for use in the PC algorithm.
define_symmetric_tests <- function(stan_model, stan_model_spike_slab) {
  bayesCItest_with_model <- function(x, y, S, suffStat, epsilon = 0.1) {
    suffStat$stan_model <- stan_model
    bayesCItest(x, y, S, suffStat, epsilon)
  }
  
  bayesCItest_spike_slab_with_model <- function(x, y, S, suffStat, epsilon = 0.1) {
    suffStat$stan_model_spike_slab <- stan_model_spike_slab
    bayesCItest_spike_slab(x, y, S, suffStat, epsilon)
  }
  
  list(
    bayesCItest_symmetric = make_symmetric_test(bayesCItest_with_model),
    bayesCItest_spike_slab_symmetric = make_symmetric_test(bayesCItest_spike_slab_with_model)
  )
}

