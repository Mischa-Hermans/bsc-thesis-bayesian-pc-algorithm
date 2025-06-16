# ---- R/stan_models.R ----
# Stan Model Compilation for Bayesian CI Testing
#
# Description:
#   Defines and compiles Stan models used for Bayesian conditional
#   independence testing in the PC algorithm. Includes both standard
#   Bayesian linear regression and spike-and-slab priors.
#
# Provides:
#   - get_stan_models(): Compiles and returns a list of Stan model objects.
#
# Dependencies:
#   - rstan
#
# Author: Mischa Hermans

get_stan_models <- function() {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Package 'rstan' must be installed.")
  }
  
  # Standard Bayesian linear regression model
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
  
  # Spike-and-slab Bayesian linear regression model
  stan_code_spike_slab <- "
    data {
      int<lower=0> N;
      int<lower=0> K;
      matrix[N, K] X;
      vector[N] y;
      vector[K] mu_beta;
      vector<lower=0>[K] sigma_beta;
      vector<lower=0, upper=1>[K] prior_inclusion;
    }
    parameters {
      real alpha;
      vector[K] beta_raw;
      vector<lower=0, upper=1>[K] inclusion_raw;
      real<lower=0> sigma;
    }
    transformed parameters {
      vector[K] beta;
      for (k in 1:K)
        beta[k] = inclusion_raw[k] * beta_raw[k];
    }
    model {
      inclusion_raw ~ beta(prior_inclusion * 10 + 1e-2, (1 - prior_inclusion) * 10 + 1e-2);
      beta_raw ~ normal(mu_beta, sigma_beta);
      y ~ normal(alpha + X * beta, sigma);
    }
  "
  
  list(
    standard = rstan::stan_model(model_code = stan_code),
    spike_slab = rstan::stan_model(model_code = stan_code_spike_slab)
  )
}
