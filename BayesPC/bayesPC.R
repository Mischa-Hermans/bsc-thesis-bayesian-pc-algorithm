# ---- bayesPC.R ----
# Bayesian PC Algorithm
#
# Description:
#   Runs the PC algorithm once using Bayesian conditional independence
#   testing via Stan-based Bayesian linear regression.
#   Allows incorporating uncertain background knowledge about causal effects:
#     - The expected effect of a variable on another can be specified in the mu matrix.
#     - The confidence in that effect can be specified via the standard deviation in the sigma matrix.
#   Outputs the estimated CPDAG and causal effects.
#
# Dependencies:
#   - rstan
#   - pcalg
#   - igraph
#
# Author: Mischa Hermans

bayesPC <- function(data,
                    alpha = 0.1,
                    mu = NULL,
                    sigma = NULL,
                    iter = 500,
                    warmup = 250,
                    adapt_delta = 0.9,
                    seed = 123) {
  
  # ---- Package Checks ----
  pkgs <- c("rstan", "pcalg", "igraph")
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' must be installed.", pkg))
    }
  }
  
  # ---- Compile Stan Model ----
  stan_code <- "
    data {
      int<lower=0> N;
      int<lower=0> K;
      matrix[N, K] X;
      vector[N] y;
      vector[K] mu_beta;
      vector<lower=0>[K] sigma_beta;
    }
    parameters {
      real alpha;
      vector[K] beta;
      real<lower=0> sigma;
    }
    model {
      beta ~ normal(mu_beta, sigma_beta);
      y ~ normal(alpha + X * beta, sigma);
    }
  "
  
  if (!exists(".__bayesPC_stan_model__", envir = .GlobalEnv)) {
    message("Compiling Stan model...")
    stan_model <- rstan::stan_model(model_code = stan_code)
    assign(".__bayesPC_stan_model__", stan_model, envir = .GlobalEnv)
  } else {
    stan_model <- get(".__bayesPC_stan_model__", envir = .GlobalEnv)
    message("Using cached Stan model.")
  }
  
  # ---- Bayesian CI Test ----
  bayesCItest <- function(x, y, S, suffStat, epsilon = 0.1) {
    data <- suffStat$data
    predictors <- if (length(S)) c(x, S) else x
    x_index <- which(predictors == x)
    
    stan_data <- list(
      N = nrow(data),
      K = length(predictors),
      X = as.matrix(data[, predictors, drop = FALSE]),
      y = data[[y]],
      mu_beta = as.array(suffStat$prior_mu[predictors, y]),
      sigma_beta = as.array(suffStat$prior_sigma[predictors, y])
    )
    
    fit <- tryCatch(
      suppressMessages(suppressWarnings(
        rstan::sampling(
          suffStat$stan_model,
          data = stan_data,
          iter = suffStat$iter,
          warmup = suffStat$warmup,
          chains = 1,
          refresh = 0,
          seed = suffStat$seed,
          control = list(adapt_delta = suffStat$adapt_delta)
        )
      )),
      error = function(e) return(NULL)
    )
    
    if (is.null(fit)) return(list(p = NA))
    
    beta_post <- rstan::extract(fit)$beta[, x_index]
    if (length(beta_post) == 0 || all(is.na(beta_post))) return(list(p = NA))
    
    list(p = 2 * min(mean(beta_post <= 0), mean(beta_post > 0)))
  }
  
  # ---- Convert Bayesian CI Test to Symmetric Form ----
  # Ensures the CI test produces the same result for X ⫫ Y | S and Y ⫫ X | S
  make_symmetric_test <- function(ci_func) {
    function(x, y, S, suffStat, alpha = 0.1) {
      cache <- suffStat$cache
      S_key <- if (length(S) == 0) "∅" else paste(sort(S), collapse = ",")
      pair_key <- paste(sort(c(x, y)), collapse = "-")
      key <- paste0(pair_key, "|", S_key)
      
      if (exists(key, envir = cache, inherits = FALSE)) {
        return(get(key, envir = cache))
      }
      
      result1 <- ci_func(x, y, S, suffStat)
      result2 <- ci_func(y, x, S, suffStat)
      p1 <- result1$p; p2 <- result2$p
      decision <- if (!is.na(p1) && !is.na(p2) && p1 > alpha && p2 > alpha) 1 else 0
      assign(key, decision, envir = cache)
      return(decision)
    }
  }
  
  bayesCItest_symmetric <- make_symmetric_test(bayesCItest)
  
  # ---- Prepare Data and Priors ----
  data <- as.data.frame(data)
  n <- nrow(data)
  p <- ncol(data)
  colnames(data) <- if (is.null(colnames(data))) paste0("X", 1:p) else colnames(data)
  
  # If no prior matrices are provided, use weakly informative priors
  if (is.null(mu)) {
    prior_mu <- matrix(0, p, p)        
  } else {
    prior_mu <- mu                      
  }
  
  if (is.null(sigma)) {
    prior_sigma <- matrix(5, p, p)     
  } else {
    prior_sigma <- sigma                
  }
  
  suffStat <- list(
    data = data,
    prior_mu = prior_mu,
    prior_sigma = prior_sigma,
    stan_model = stan_model,
    iter = iter,
    warmup = warmup,
    adapt_delta = adapt_delta,
    seed = seed,
    cache = new.env(hash = TRUE, parent = emptyenv())
  )
  
  # ---- Run PC Algorithm ----
  message("Running Bayesian PC algorithm...")
  pc_fit <- pcalg::pc(
    suffStat = suffStat,
    indepTest = bayesCItest_symmetric,
    alpha = alpha,
    labels = colnames(data)
  )
  
  g_est <- igraph::graph_from_graphnel(pc_fit@graph)
  adj_est <- as.matrix(igraph::as_adjacency_matrix(g_est, sparse = FALSE))
  
  # ---- Estimate Effects (OLS per parent set) ----
  estimate_effects <- function(data, adj_matrix) {
    p <- ncol(data)
    effects <- matrix(NA, p, p)
    colnames(effects) <- rownames(effects) <- colnames(data)
    
    for (j in 1:p) {
      parents <- which(adj_matrix[, j] == 1)
      if (length(parents) == 0) next
      
      y <- data[[j]]
      X_parents <- as.data.frame(data[, parents, drop = FALSE])
      
      # Give predictor columns proper names
      colnames(X_parents) <- colnames(data)[parents]
      
      # Fit regression y ~ parent predictors
      fit <- lm(y ~ ., data = X_parents)
      
      # Store coefficients (excluding intercept)
      coefs <- coef(fit)[-1]
      effects[parents, j] <- coefs
    }
    
    effects
  }
  
  effects <- estimate_effects(data, adj_est)
  
  message("Bayesian PC completed.")
  list(
    graph = g_est,
    adjacency = adj_est,
    effects = effects
  )
}
