library(microbenchmark)

# --- Load both implementations ----------------------------------------------
# R reference implementation and its helpers
devtools::load_all()

# --- Data -------------------------------------------------------------------
# Synthetic latent structure as in tests/testthat helpers
K <- 4
R <- 4
npc <- 30
n1 <- npc * K # number of rows
n2 <- 30 # number of columns
sigma2 <- 1
niter_metropolis <- 1L # single Metropolis step per Gibbs iteration

sigma2_fixed <- TRUE
minibatch <- TRUE

rho_within <- 0.9
rho_between <- -1 / (K * npc - 1)

set.seed(1)
Sigma <- matrix(rho_between, nrow = n1, ncol = n1)
same_cluster <- outer(rep(seq_len(K), each = npc), rep(seq_len(K), each = npc), `==`)
Sigma[same_cluster] <- rho_within
diag(Sigma) <- 1

#
set.seed(2)
P <- rmatrixnormal(M = matrix(0, n1, K - 1), U = Sigma, V = diag(1, K - 1))
set.seed(2)
Pcpp <- rmatrixnormal_cpp(M = matrix(0, n1, K - 1), U = Sigma, V = diag(1, K - 1))
P - Pcpp

# --- Benchmark 1: standalone P samplers (one Metropolis pass) ---------------
bm_rmatrixnormal <- microbenchmark(
  "R" = rmatrixnormal(M = matrix(0, n1, K - 1), U = Sigma, V = diag(1, K - 1)),
  "C++" = rmatrixnormal_cpp(M = matrix(0, n1, K - 1), U = Sigma, V = diag(1, K - 1)),
  times = 100,
  unit = "ms"
)

bm_rmatrixnormal_summary <- summary(bm_rmatrixnormal)
cat("\n--- Benchmark 1: sampling of Matrix Normal ---\n")
print(bm_rmatrixnormal_summary)

median_rmatrixnormal <- bm_rmatrixnormal_summary$median
speedup <- median_rmatrixnormal[1] / median_rmatrixnormal[2]
cat(sprintf(
  "Median: R = %.2f ms | C++ = %.2f ms | C++ classical is %.4fx faster\n",
  median_rmatrixnormal[1], median_rmatrixnormal[2], speedup
))


# --- Benchmark 2: density computation ---------------
bm_dmatrixnormal <- microbenchmark(
  "R plain" = dmatrixnormal_plain(P, M = matrix(0, n1, K - 1), U = Sigma, V = diag(1, K - 1)),
  "R Schur" = dmatrixnormal_schur(P, M = matrix(0, n1, K - 1), U = Sigma, V = diag(1, K - 1)),
  "C++" = dmatrixnormal_cpp(P, M = matrix(0, n1, K - 1), U = Sigma, V = diag(1, K - 1)),
  times = 1000,
  unit = "ms"
)

bm_dmatrixnormal_summary <- summary(bm_dmatrixnormal)
cat("\n--- Benchmark 2: density of Matrix Normal ---\n")
print(bm_dmatrixnormal_summary)

median_dmatrixnormal <- bm_dmatrixnormal_summary$median
speedup <- median_dmatrixnormal[1] / median_dmatrixnormal[2]
cat(sprintf(
  "Median: R = %.2f ms | C++ = %.2f ms | C++ classical is %.3fx faster\n",
  median_dmatrixnormal[1], median_dmatrixnormal[2], speedup
))

# --- Plots (optional, requires ggplot2) --------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  bm_all <- rbind(
    transform(as.data.frame(bm_mh), level = "P sampler only"),
    transform(as.data.frame(bm_gibbs), level = "Gibbs chain")
  )
  bm_all$expr <- factor(bm_all$expr, levels = names(samplers))
  g <- ggplot2::ggplot(bm_all, ggplot2::aes(x = expr, y = time / 1e6, fill = expr)) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_wrap(~level, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Time (ms)", title = "C++ vs R Metropolis-Hastings P samplers") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 15, hjust = 1))
  gg_path <- file.path("dev", "bench_cpp_metropolis_gibbs.png")
  ggplot2::ggsave(gg_path, g, width = 8, height = 4)
  cat("Boxplot saved to:", gg_path, "\n")
}
