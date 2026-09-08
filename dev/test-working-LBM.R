devtools::load_all()
library(sbm)
library(posterior)
library(bayesplot)
library(future)
plan("multicore")

palette_okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
color_scheme_set(scheme = palette_okabe_ito[1:6])


eps <- 0.05
true_pi <- c(0.1, 0.3, 0.6) # c(0.25, 0.3, 0.45)
true_rho <- c(0.1, 0.3, 0.6)
nbNodes <- c(100, 100)
K <- length(true_pi)
R <- length(true_rho)

true_alpha <- matrix(seq(1, 10 * K * R, by = 10), nrow = K)

bipartite_sampler <- sampleBipartiteSBM(nbNodes = nbNodes, blockProp = list(true_pi, true_rho), connectParam = list(mean = true_alpha), model = "poisson")
set.seed(123)
lbmRealisation <- bipartite_sampler$rNetwork()
Y <- lbmRealisation[["networkData"]]
true_Z_ind <- lbmRealisation[["indMemberships"]][[1]]
true_W_ind <- lbmRealisation[["indMemberships"]][[2]]

check_lbm_identifiability(netMat = Y, alpha = true_alpha, pi = true_pi, rho = true_rho, K = K, R = R)

priors_hyper_params <- list(
    etas_0 = rep(2, K),
    gammas_0 = rep(2, R),
    a0 = 1, b0 = 1 / 30
)

alpha_recovery_graphs <- function(...) {
    alpha_recovery_plots(..., true_alpha = true_alpha, find_permutation = find_permutation_alphas_L2)
}

rho_recovery_graphs <- function(...) {
    rho_recovery_plots(..., true_alpha = true_alpha, true_rho = true_rho, true_W_ind = true_W_ind)
}

mean_pi_recovery_graphs <- function(...) {
    mean_pi_recovery_plots(..., true_alpha = true_alpha, true_Z_ind = true_Z_ind)
}

pi_recovery_graphs <- function(...) {
    pi_recovery_plots(..., true_alpha = true_alpha, true_Z_ind = true_Z_ind)
}



nbiter <- 500L
nbchains <- 4
## Classical LBM

lbm_results <- chains_gibbs_sampling_lbm_poisson(nchains = nbchains, Y = Y, niter = nbiter, priors_hyper_params = priors_hyper_params, init_Z = NULL, init_W = NULL, K, R)

lbm_stan <- lbm_results_to_stan_draws(lbm_results, K = K, R = R, apply_burnin_thinning = FALSE)


pi_recovery_graphs(lbm_stan) + patchwork::plot_annotation(title = "Classical LBM only")
rho_recovery_graphs(lbm_stan) + patchwork::plot_annotation(title = "Classical LBM only")
alpha_recovery_graphs(lbm_stan) + patchwork::plot_annotation(title = "Classical LBM only")

# no_indicator_variables <- grepl(variables(lbm_stan), pattern = "^[^ZW][a-zA-Z]*\\[[0-9]*")

# summarise_draws(lbm_stan[, , no_indicator_variables]) |> View()
# summarise_draws(subset_draws(lbm_stan, chains = 1)[, , no_indicator_variables])


ARI_table(lbm_stan, true = true_Z_ind) |> colMeans()
ARI_table(lbm_stan, true = true_W_ind, "W") |> colMeans()


mcmc_trace(lbm_stan, regex_pars = "^[^ZW][a-zA-Z]*.[0-9]*") + patchwork::plot_annotation(title = "Classical LBM only")

## Classical LBM forced_order

lbm_forced_results <- chains_gibbs_sampling_lbm_poisson(nchains = nbchains, Y = Y, niter = nbiter, priors_hyper_params = priors_hyper_params, init_Z = NULL, init_W = NULL, K, R, force_order = TRUE)

lbm_forced_stan <- lbm_results_to_stan_draws(lbm_forced_results, K = K, R = R, apply_burnin_thinning = FALSE)
pi_recovery_graphs(lbm_forced_stan) + patchwork::plot_annotation(title = "Classical LBM forced order")
rho_recovery_graphs(lbm_forced_stan) + patchwork::plot_annotation(title = "Classical LBM forced order")
alpha_recovery_graphs(lbm_forced_stan) + patchwork::plot_annotation(title = "Classical LBM forced order")

ARI_table(lbm_forced_stan, true = true_Z_ind) |> colMeans()
ARI_table(lbm_forced_stan, true = true_W_ind, "W") |> colMeans()


mcmc_trace(lbm_forced_stan, regex_pars = "^[^ZW][a-zA-Z]*.[0-9]*") + patchwork::plot_annotation(title = "Classical LBM forced order")

##  Classical LBM with sbm init

fitted_lbm <- estimateBipartiteSBM(netMat = Y, model = "poisson")
fitted_lbm_Z <- onehot_encode(fitted_lbm[["memberships"]][[1]], K = K)
fitted_lbm_W <- onehot_encode(fitted_lbm[["memberships"]][[2]], K = R)
lbm_with_init_results <- chains_gibbs_sampling_lbm_poisson(nchains = nbchains, Y = Y, niter = nbiter, priors_hyper_params = priors_hyper_params, init_Z = fitted_lbm_Z, init_W = fitted_lbm_W, K, R)

lbm_with_init_stan <- lbm_results_to_stan_draws(lbm_with_init_results, K = K, R = R)

pi_recovery_graphs(lbm_with_init_stan) + patchwork::plot_annotation(title = "Classical LBM with SBM init")
rho_recovery_graphs(lbm_with_init_stan) + patchwork::plot_annotation(title = "Classical LBM with SBM init")
alpha_recovery_graphs(lbm_with_init_stan) + patchwork::plot_annotation(title = "Classical LBM with SBM init")

ARI_table(lbm_with_init_stan, true = true_Z_ind) |> colMeans()
ARI_table(lbm_with_init_stan, true = true_W_ind, "W") |> colMeans()


mcmc_trace(lbm_with_init_stan, regex_pars = "^[^ZW][a-zA-Z]*.[0-9]*") + patchwork::plot_annotation(title = "Classical LBM with SBM init")

##  LBM with an Identity latent cov

### sigma2 = 1
lbmcovs1_results <- chains_gibbs_sampling_lbm_cov_poisson(nchains = nbchains, Y = Y, niter = nbiter, niter_metropolis = 1, sigma2_fixed = 1, Sigma = diag(rep(1, nbNodes[1])), init_Z = NULL, init_W = NULL, K = K, R = R, tol = TOL, P_sampler = sample_P_metropolis_classical_cpp, priors_hyper_params = priors_hyper_params)



lbmcovs1_stan <- lbm_results_to_stan_draws(lbmcovs1_results, K = K, R = R, apply_burnin_thinning = TRUE, delabel_switch = FALSE)


lbmcovs1_stan <- delabel_switch_stan_true(draws = lbmcovs1_stan, K, R, true_alpha = true_alpha)

mean_pi_recovery_graphs(lbmcovs1_stan) + patchwork::plot_annotation(title = "LBM covar only sigma2=1")
rho_recovery_graphs(lbmcovs1_stan) + patchwork::plot_annotation(title = "LBM covar only sigma2=1")
alpha_recovery_graphs(lbmcovs1_stan) + patchwork::plot_annotation(title = "LBM covar only sigma2=1")

ARI_table(lbmcovs1_stan, true = true_Z_ind)
ARI_table(lbmcovs1_stan, true = true_W_ind, "W") |> colMeans()

draws <- lbmcovs1_stan
draws_Z <- subset_draws(draws, variable = "Z", iteration = 500, chain = 3) |>
    as.vector() |>
    .one_hot(K)

draws_W <- subset_draws(draws, variable = "W", iteration = 500, chain = 3) |>
    as.vector() |>
    .one_hot(R)

A0 <- t(draws_Z) %*% Y %*% draws_W
B0 <- t(draws_Z) %*% matrix(1, nbNodes[1], nbNodes[2]) %*% draws_W


(1 + A0) %/% (B0 + 1 / 30)


mcmc_trace(lbmcovs1_stan, regex_pars = "^[^PZW][a-zA-Z]*.[0-9]*") + patchwork::plot_annotation(title = "LBM covar only sigma2=1")

### sigma2 = 5
lbmcovs5_results <- chains_gibbs_sampling_lbm_cov_poisson(nchains = nbchains, Y = Y, niter = nbiter, niter_metropolis = 1, sigma2_fixed = 5, Sigma = diag(rep(1, nbNodes[1])), init_Z = NULL, init_W = NULL, K = K, R = R, tol = TOL, P_sampler = sample_P_metropolis_classical_cpp)


lbmcovs5_stan <- lbm_results_to_stan_draws(lbmcovs5_results, K = K, R = R)

mean_pi_recovery_graphs(lbmcovs5_stan) + patchwork::plot_annotation(title = "LBM covar only sigma2=5")
rho_recovery_graphs(lbmcovs5_stan) + patchwork::plot_annotation(title = "LBM covar only sigma2=5")
alpha_recovery_graphs(lbmcovs5_stan) + patchwork::plot_annotation(title = "LBM covar only sigma2=5")

ARI_table(lbmcovs5_stan, true = true_Z_ind) |> colMeans()
ARI_table(lbmcovs5_stan, true = true_W_ind, "W") |> colMeans()

### sigma2 = 200
lbmcovs200_results <- chains_gibbs_sampling_lbm_cov_poisson(nchains = nbchains, Y = Y, niter = nbiter, niter_metropolis = 1, sigma2_fixed = 200, Sigma = diag(rep(1, nbNodes[1])), init_Z = NULL, init_W = NULL, K = K, R = R, tol = TOL, P_sampler = sample_P_metropolis_classical_cpp)


lbmcovs200_stan <- lbm_results_to_stan_draws(lbmcovs200_results, K = K, R = R)

mean_pi_recovery_graphs(lbmcovs200_stan) + patchwork::plot_annotation(title = "LBM covar only sigma2=200")
rho_recovery_graphs(lbmcovs200_stan) + patchwork::plot_annotation(title = "LBM covar only sigma2=200")
alpha_recovery_graphs(lbmcovs200_stan) + patchwork::plot_annotation(title = "LBM covar only sigma2=200")

ARI_table(lbmcovs200_stan, true = true_Z_ind) |> colMeans()
ARI_table(lbmcovs200_stan, true = true_W_ind, "W") |> colMeans()

### sigma2 = 500
lbmcovs500_results <- chains_gibbs_sampling_lbm_cov_poisson(nchains = nbchains, Y = Y, niter = nbiter, niter_metropolis = 1, sigma2_fixed = 500, Sigma = diag(rep(1, nbNodes[1])), init_Z = NULL, init_W = NULL, K = K, R = R, tol = TOL, P_sampler = sample_P_metropolis_classical_cpp)


lbmcovs500_stan <- lbm_results_to_stan_draws(lbmcovs500_results, K = K, R = R)

mean_pi_recovery_graphs(lbmcovs500_stan) + patchwork::plot_annotation(title = "LBM covar only sigma2=500")
rho_recovery_graphs(lbmcovs500_stan) + patchwork::plot_annotation(title = "LBM covar only sigma2=500")
alpha_recovery_graphs(lbmcovs500_stan) + patchwork::plot_annotation(title = "LBM covar only sigma2=500")

ARI_table(lbmcovs500_stan, true = true_Z_ind) |> colMeans()
ARI_table(lbmcovs500_stan, true = true_W_ind, "W") |> colMeans()

##  LBM with an Identity latent cov and spectral Z init
spec_Z <- .one_hot(colSBM:::spectral_clustering(X = Y, K = K), K)

lbmcov_with_spectralZinit_results <- chains_gibbs_sampling_lbm_cov_poisson(nchains = nbchains, Y = Y, niter = nbiter, niter_metropolis = 1, sigma2_fixed = TRUE, Sigma = diag(rep(1, nbNodes[1])), init_Z = spec_Z, init_W = NULL, K = K, R = R, tol = TOL, P_sampler = sample_P_metropolis_classical_cpp)


lbmcov_with_spectralZinit_stan <- lbm_results_to_stan_draws(lbmcov_with_spectralZinit_results, K = K, R = R)

mean_pi_recovery_graphs(lbmcov_with_spectralZinit_stan) + patchwork::plot_annotation(title = "LBM covar with spectral Z init")
rho_recovery_graphs(lbmcov_with_spectralZinit_stan) + patchwork::plot_annotation(title = "LBM covar with spectral Z init")
alpha_recovery_graphs(lbmcov_with_spectralZinit_stan) + patchwork::plot_annotation(title = "LBM covar with spectral Z init")

mcmc_trace(lbmcov_with_spectralZinit_stan, regex_pars = "^[^ZWP][a-zA-Z]*.[0-9]*") + patchwork::plot_annotation(title = "LBM covar with spectral Z init")


ARI_table(lbmcov_with_spectralZinit_stan, true = true_Z_ind) |> colMeans()
ARI_table(lbmcov_with_spectralZinit_stan, true = true_W_ind, "W") |> colMeans()

##  LBM with an Identity latent cov and spectral init
spec_Z <- .one_hot(colSBM:::spectral_clustering(X = Y, K = K), K)
spec_W <- .one_hot(colSBM:::spectral_clustering(X = t(Y), K = R), R)

lbmcov_with_spectralinit_results <- chains_gibbs_sampling_lbm_cov_poisson(nchains = nbchains, Y = Y, niter = nbiter, niter_metropolis = 1, sigma2_fixed = TRUE, Sigma = diag(rep(1, nbNodes[1])), init_Z = spec_Z, init_W = spec_W, K = K, R = R, tol = TOL, P_sampler = sample_P_metropolis_classical_cpp)


lbmcov_with_spectralinit_stan <- lbm_results_to_stan_draws(lbmcov_with_spectralinit_results, K = K, R = R)

mean_pi_recovery_graphs(lbmcov_with_spectralinit_stan) + patchwork::plot_annotation(title = "LBM covar with spectral init")
rho_recovery_graphs(lbmcov_with_spectralinit_stan) + patchwork::plot_annotation(title = "LBM covar with spectral init")
alpha_recovery_graphs(lbmcov_with_spectralinit_stan) + patchwork::plot_annotation(title = "LBM covar with spectral init")

mcmc_trace(lbmcov_with_spectralinit_stan, regex_pars = "^[^ZWP][a-zA-Z]*.[0-9]*") + patchwork::plot_annotation(title = "LBM covar with spectral init")


ARI_table(lbmcov_with_spectralZinit_stan, true = true_Z_ind) |> colMeans()
ARI_table(lbmcov_with_spectralZinit_stan, true = true_W_ind, "W") |> colMeans()

##  LBM with an Identity latent cov and SBM Z init
lbmcov_with_sbmZinit_results <- chains_gibbs_sampling_lbm_cov_poisson(nchains = nbchains, Y = Y, niter = nbiter, niter_metropolis = 1, sigma2_fixed = 50, Sigma = diag(rep(1, nbNodes[1])), init_Z = fitted_lbm_Z, init_W = NULL, K = K, R = R, tol = TOL, P_sampler = sample_P_metropolis_classical_cpp)


lbmcov_with_sbmZinit_stan <- lbm_results_to_stan_draws(lbmcov_with_sbmZinit_results, K = K, R = R)

mean_pi_recovery_graphs(lbmcov_with_sbmZinit_stan) + patchwork::plot_annotation(title = "LBM covar with SBM Z init")
rho_recovery_graphs(lbmcov_with_sbmZinit_stan) + patchwork::plot_annotation(title = "LBM covar with SBM Z init")
alpha_recovery_graphs(lbmcov_with_sbmZinit_stan) + patchwork::plot_annotation(title = "LBM covar with SBM Z init")

mcmc_trace(lbmcov_with_sbmZinit_stan, regex_pars = "^[^ZWP][a-zA-Z]*.[0-9]*") + patchwork::plot_annotation(title = "LBM covar with SBM init")

ARI_table(lbmcov_with_sbmZinit_stan, true = true_Z_ind) |> colMeans()
ARI_table(lbmcov_with_sbmZinit_stan, true = true_W_ind, "W") |> colMeans()


##  LBM with an Identity latent cov and SBM init
lbmcov_with_sbminit_results <- chains_gibbs_sampling_lbm_cov_poisson(nchains = nbchains, Y = Y, niter = nbiter, niter_metropolis = 1, sigma2_fixed = TRUE, Sigma = diag(rep(1, nbNodes[1])), init_Z = fitted_lbm_Z, init_W = fitted_lbm_W, K = K, R = R, tol = TOL, P_sampler = sample_P_metropolis_classical_cpp, priors_hyper_params = priors_hyper_params)


lbmcov_with_sbminit_stan <- lbm_results_to_stan_draws(lbmcov_with_sbminit_results, K = K, R = R)

mean_pi_recovery_graphs(lbmcov_with_sbminit_stan) + patchwork::plot_annotation(title = "LBM covar with SBM init")
rho_recovery_graphs(lbmcov_with_sbminit_stan) + patchwork::plot_annotation(title = "LBM covar with SBM init")
alpha_recovery_graphs(lbmcov_with_sbminit_stan) + patchwork::plot_annotation(title = "LBM covar with SBM init")

mcmc_trace(lbmcov_with_sbminit_stan, regex_pars = "^[^ZWP][a-zA-Z]*.[0-9]*") + patchwork::plot_annotation(title = "LBM covar with SBM init")
