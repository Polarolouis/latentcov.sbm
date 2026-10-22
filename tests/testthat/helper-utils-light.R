niter <- 100
n1 <- 30
n2 <- 50

K <- 3
R <- 4
set.seed(42)
Z_draws <- sapply(seq(n1), function(i) {
  sample.int(n = K, size = niter, replace = TRUE, prob = c(0.1, 0.7, 0.2))
})
colnames(Z_draws) <- paste0("Z[", seq(n1), "]")
Z_draws <- posterior::as_draws_array(Z_draws)
Z_posterior_probs <-
  apply(Z_draws, 3, function(x) {
    prop.table(table(x))
  }) |>
  t()

true_Z <- factor(apply(Z_posterior_probs, 1, which.max), levels = seq(K))
W_draws <- sapply(seq(n2), function(j) {
  sample.int(n = R, size = niter, replace = TRUE, prob = c(0.1, 0.2, 0.5, 0.1))
})
colnames(W_draws) <- paste0("W[", seq(n2), "]")
W_draws <- posterior::as_draws_array(W_draws)
W_posterior_probs <- apply(W_draws, 3, function(x) {
  prop.table(table(x))
}) |>
  t()

true_W <- factor(apply(W_posterior_probs, 1, which.max), levels = seq(R))

draws <- posterior::bind_draws(Z_draws, W_draws)
