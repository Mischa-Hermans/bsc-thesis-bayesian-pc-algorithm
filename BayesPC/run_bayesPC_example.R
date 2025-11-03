source("bayesPC.R")

set.seed(123)

# Generate some data
n <- 200
X1 <- rnorm(n)
X2 <- 2 * X1 + rnorm(n)
X3 <- -1.5 * X1 + rnorm(n)
X4 <- 0.5 * X2 - 0.8 * X3 + rnorm(n)
X5 <- 3 * X1 + rnorm(n)
data <- data.frame(X1, X2, X3, X4, X5)

p <- ncol(data)

# ---- Define prior expected effects ----
mu <- matrix(0, p, p) # No effect if no prior specified
mu[1, 2] <- 2      # X1 -> X2
mu[1, 3] <- -1.5   # X1 -> X3
mu[1, 5] <- 3      # X1 -> X5
mu[2, 4] <- 0.5    # X2 -> X4
mu[3, 4] <- -0.8   # X3 -> X4

# ---- Define confidence in priors ----
sigma <- matrix(5, p, p) # Large variance if no prior specified    
sigma[1, 2] <- 0.1
sigma[1, 3] <- 0.1   
sigma[1, 5] <- 0.1      
sigma[2, 4] <- 0.1
sigma[3, 4] <- 0.1

# ---- Run Bayesian PC ----
result <- bayesPC(data, alpha = 0.2, mu = mu, sigma = sigma)

# ---- Show results ----
cat("Estimated causal effects:\n")
print(result$effects)

cat("\nEstimated adjacency matrix:\n")
print(result$adjacency)

cat("Plotting estimated causal graph...\n")
plot(result$graph, vertex.label = colnames(data), edge.arrow.size = 0.5)
