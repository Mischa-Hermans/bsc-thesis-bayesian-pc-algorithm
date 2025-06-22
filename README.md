# Thesis Code: Incorporating Background Knowledge in Graphical Causal Models

This repository contains all simulation code for the bachelor's thesis:

**Title**: *Incorporating Background Knowledge in Graphical Causal Models: A Bayesian Approach*  
**Author**: Mischa Hermans  
**Institution**: Maastricht University – BSc Econometrics & Operations Research  
**Supervisor**: Dr. Nalan Bastürk

---

## Overview

This project compares different ways of integrating prior knowledge into the PC algorithm for causal structure learning:

- **No Background**: standard PC algorithm
- **Hard Constraints**: background knowledge enforced as fixed edges/gaps
- **Bayesian Priors**: CI tests informed by priors on regression coefficients using Stan
- **Spike-and-Slab**: CI tests using sparse Bayesian priors on inclusion

---

## Repository Structure

```
bsc-thesis-bayesian-pc-algorithm/
├── models/                     # Stan model compilation
│   └── stan_models.R
│
├── sensitivity/                # Sensitivity analysis framework
│   ├── main_sensitivity_analysis.R
│   └── sensitivity_analysis_functions.R
│
├── simulation/                 # Simulation runner and logic
│   ├── main_simulation.R
│   └── simulation_functions.R
│
├── utils/                      # Utility functions
│   ├── ci_tests.R              # Conditional independence tests
│   ├── data_generation.R       # Data simulation from DAG
│   ├── evaluation_metrics.R    # Confusion matrix and effect metrics
│   └── prior_generation.R      # Sampling prior edges and gaps
│
├── .gitignore                  # Git ignore file
├── LICENSE                     # MIT License
└── README.md                   # Project description and instructions

```

---

## Requirements

Install required R packages with:

```r
install.packages(c(
  "pcalg", "igraph", "rstan", "ggplot2", "dplyr", "tidyr",
  "openxlsx", "kableExtra", "purrr"
))
```

---

## Running the Code

From the project root, simply run:

```r
source("simulation/main_simulation.R")
```

This script will:
- Generate simulated data from a fixed DAG
- Run the PC algorithm with different prior configurations
- Evaluate structure recovery and causal effect accuracy by calculating Confusion metrics and MSE and variance of estimated causal effects

To run the sensitivity analysis, use:

```r
source("sensitivity/main_sensitivity_analysis.R")
```

This script will:

- Run multiple simulation rounds while varying one parameter at a time (e.g., pTP, pFP, mu, sigma, n, or noise_sd)
- Evaluate how structure recovery and causal effect estimation performance change with each parameter
- Print summary tables for confusion metrics and MSE/variance of estimated effects across all settings and parameter values

---

---

---

## Acknowledgements

The Bayesian CI test framework builds on Basturk et al. (2024), extended with Stan-based regression priors and spike-and-slab modeling.

=======
# bsc-thesis-causal-discovery-with-priors
Bayesian extensions of the PC algorithm for causal discovery using prior knowledge. Bachelor thesis project in Econometrics &amp; Operations Research at Maastricht University.