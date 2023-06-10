source("R/average_partial_effect.R")
source("R/minimax.R")
source("R/plugin.R")
source("R/rlasso.R")
# source("R/minimax.R")
set.seed(1)

n = 300; p = 600; nclust = 10
beta = 2 / (1:p) / sqrt(sum(1/(1:p)^2))
clust.alpha = rep(c(0.3, 0.7), nclust/2)
cluster.center = 0.5 * matrix(rnorm(nclust * p), nclust, p)
cluster = sample.int(nclust, n, replace = TRUE)

X = cluster.center[cluster,] + matrix(rnorm(n * p), n, p)
W = rbeta(n, clust.alpha[cluster], 1 - clust.alpha[cluster])
Y = X %*% beta + rnorm(n, 0, 1) + 2 * W * X[,2]

tau.hat = average_partial_effect(X, Y, W, solver="MOSEK")
print(paste("true tau:", round(mean(2 * cluster.center[,2]), 2)))
print(paste("point estimate:", round(tau.hat[1], 2)))
print(paste0("95% CI for tau: (", round(tau.hat[1] - 1.96 * tau.hat[2], 2), ", ", round(tau.hat[1] + 1.96 * tau.hat[2], 2), ")"))
