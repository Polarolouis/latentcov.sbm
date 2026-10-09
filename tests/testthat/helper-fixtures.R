# Shared fixtures for unit tests.
# Small, deterministic objects reused across test files.

K_fx <- 3L
R_fx <- 2L
n1_fx <- 6L
n2_fx <- 5L

# Small block-structured correlation matrix (6 x 6), 2 blocks of 3.
block_members <- rep(1:2, each = 3)
Sigma_fx <- matrix(-0.3, n1_fx, n1_fx)
Sigma_fx[outer(block_members, block_members, `==`)] <- 0.6
diag(Sigma_fx) <- 1

# Deterministic small latent positions (6 x K-1 = 6 x 2).
P_fx <- matrix(c(
  1.2, -0.4,
  -0.6, 0.9,
  0.1, 0.3,
  -0.9, -1.1,
  0.7, 0.2,
  -0.3, 1.1
), nrow = n1_fx, byrow = TRUE)
colnames(P_fx) <- paste0("coord", seq_len(K_fx - 1L))

# One-hot row memberships (n1 x K).
Z_fx <- matrix(0, n1_fx, K_fx)
Z_fx[cbind(seq_len(n1_fx), c(1L, 2L, 3L, 1L, 2L, 3L))] <- 1

# One-hot column memberships (n2 x R).
W_fx <- matrix(0, n2_fx, R_fx)
W_fx[cbind(seq_len(n2_fx), c(1L, 2L, 1L, 2L, 1L))] <- 1

# Connectivity matrix.
alpha_fx <- matrix(c(3, 1, 2, 1, 2, 4), nrow = K_fx, ncol = R_fx, byrow = TRUE)

# Row / column block proportions.
pi_fx <- c(0.5, 0.3, 0.2)
rho_fx <- c(0.4, 0.6)

# Deterministic Poisson mesh (n1 x n2) with a couple of missings.
Y_fx <- matrix(c(
  2, 3, 1, 2, 0,
  1, NA, 3, 1, 1,
  4, 2, 1, 5, 6,
  NA, 1, 0, 1, 2,
  1, 4, 3, 3, 1,
  0, 1, 2, 2, 4
), nrow = n1_fx, byrow = TRUE)
mask_fx <- (!is.na(Y_fx)) * 1L

# A tiny positive-definite Sigma for sampler-level tests.
Sigma_spd_fx <- diag(6) * 1 + 0.1
