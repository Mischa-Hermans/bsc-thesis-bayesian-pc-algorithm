<<<<<<< HEAD
# 🧠 Thesis Code: Incorporating Background Knowledge in Graphical Causal Models

This repository contains all simulation code for the bachelor's thesis:

**Title**: *Incorporating Background Knowledge in Graphical Causal Models: A Bayesian Approach*  
**Author**: Mischa Hermans 
**Institution**: Maastricht University – BSc Econometrics & Operations Research  
**Supervisor**: Dr. Nalan Basturk

---

## 📘 Overview

This project compares different ways of integrating prior knowledge into the PC algorithm for causal structure learning:

- **No Background**: standard PC algorithm
- **Hard Constraints**: background knowledge enforced as fixed edges/gaps
- **Bayesian Priors**: CI tests informed by priors on regression coefficients using Stan
- **Spike-and-Slab**: CI tests using sparse Bayesian priors on inclusion
- **Add Data** *(optional)*: injects fake data to reinforce prior structure

---

## 📁 Repository Structure

```
thesis-causal-discovery/
├── R/                        # Modular R code
│   ├── stan_models.R         # Stan model compilation
│   ├── ci_tests.R            # CI testing functions
│   ├── simulation_helpers.R  # Prior sampling, dummy data
│   ├── evaluation.R          # Confusion and effect tracking
│   └── plotting.R            # Sensitivity plots
├── models/                   # (Optional) .stan model files
│   ├── ci_test.stan
│   └── ci_spike_slab.stan
├── results/                  # Output Excel files and plots
├── run_simulation.R          # Main simulation script
└── README.md                 # This file
```

---

## 📦 Requirements

Install required R packages with:

```r
install.packages(c(
  "pcalg", "igraph", "rstan", "ggplot2", "dplyr", "tidyr",
  "openxlsx", "kableExtra", "purrr"
))
```

Also ensure `rstan` is correctly configured to compile Stan models.

---

## ▶️ Running the Code

From the project root, simply run:

```r
source("run_simulation.R")
```

This script will:
- Generate simulated data from a fixed DAG
- Run the PC algorithm with different prior configurations
- Evaluate structure recovery and causal effect accuracy
- Save results to `results/simulation_results.xlsx`
- Plot confusion and effect metrics

---

## 📊 Output

You’ll get:
- Confusion metrics (TP, FP, FN, TN) by variable category
- MSE and variance of estimated causal effects
- Sensitivity plots saved to PDF (if desired)

---

## 📜 License

This project is licensed under the MIT License.

---

## 🧠 Acknowledgements

The Bayesian CI test framework builds on Basturk et al. (2021), extended with Stan-based regression priors and spike-and-slab modeling.

=======
# bsc-thesis-causal-discovery-with-priors
Bayesian extensions of the PC algorithm for causal discovery using prior knowledge. Bachelor thesis project in Econometrics &amp; Operations Research at Maastricht University.
>>>>>>> ff49f1dc0ae285b4e481fbc514aac9f381efa044
