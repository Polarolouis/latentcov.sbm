# Common
seed <- 1234L
set.seed(seed)
# Rows
K <- 3
sigma2star <- 1
npc <- 30
n1 <- npc * K
phylo_memberships <- rep(seq(K), each = npc)
rho_within <- 0.9
rho_between <- -1 / (K * npc - 1)

# indep
Pindep <- matrix(rnorm(n1 * (K - 1)), nrow = n1)

# non indep
library(mvtnorm)

Sigma <- matrix(rho_between, nrow = n1, ncol = n1)
same_cluster <- outer(phylo_memberships, phylo_memberships, `==`)
Sigma[same_cluster] <- rho_within
diag(Sigma) <- 1


Pnonindep <- simulate_P_and_Z(K = K, Sigma, sigma2star, 1)[["P"]]
