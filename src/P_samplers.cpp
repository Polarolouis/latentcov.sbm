#include "P_samplers.h"
#include "mvnorm.h"
#include "utils.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// Uniform-random permutation drawn from R's RNG (not Armadillo's). Mirrors
// arma::randperm(n): draw `n` uniforms and argsort them, so the result and the
// RNG consumption match `order(runif(n))` in R.
inline arma::uvec random_perm_r(int n) {
  arma::vec u = Rcpp::runif(n);
  return arma::sort_index(u);
}
//' A Metropolis-Hastings sampler for P with a classical random walk (C++)
//'
//' Proposes each row \eqn{P_i} with an isotropic Gaussian random walk of
//' variance \code{rho} and accepts with the full Metropolis-Hastings ratio
//' combining the multinomial prior and the \eqn{N(\mu_i, \Sigma_i)}
//' conditional.
//' @param P a matrix of size \eqn{n_1 \times K-1} specifying latent
//'   position of the nodes probabilities of membership in the latent space
//' @param Z a matrix of size \eqn{n_1\times K} with a single 1 per line
//'   indicating the membership of row node \eqn{i}, \eqn{Z_{i,k} = 1} if
//'   \eqn{i} is in group \eqn{k} 0 else
//' @param Sigma a covariance matrix describing the covariance between row
//'   nodes \eqn{(n_1\times n_1)}
//' @param sigma2 a variance parameter indicating the variance between the
//'   K-1 columns of P
//' @param minibatch a boolean indicating wether to update the rows in the
//'   lexical order or to sample at each iteration. Default to TRUE
//' @param niter_metropolis an integer specifying the number of metropolis
//'   iterations to perform. Defaults to 50.
//' @param rho a double indicating the variance of the gaussian proposal
//'   distribution.
//' @return an updated matrix of latent positions, of the same size as
//'   \code{P}
//' @seealso [sample_P_metropolis_classical()],
//' [sample_P_metropolis_trick()],
//' [sample_P_metropolis_trick_cpp()]
//' @export
// [[Rcpp::export]]
arma::mat sample_P_metropolis_classical_cpp(arma::mat &P, arma::mat &Z,
                                            arma::mat &Sigma, double sigma2,
                                            bool minibatch = true,
                                            int niter_metropolis = 50,
                                            double rho = 1.0) {
  Rcpp::RNGScope scope;
  mat running_P = P;
  int n = running_P.n_rows;
  int K = running_P.n_cols + 1;
  arma::mat basis = default_Psi_function_cpp(K);

  arma::rowvec rho_values = {1.0, 0.1, 10};

  arma::uvec row_order;

  if (minibatch) {
    row_order = arma::randperm(n);
  } else {
    row_order = arma::regspace<uvec>(0, n - 1);
  }

  // Loop over the individuals
  int accepted_count = 0;
  arma::uvec::iterator ind_end = row_order.end();
  for (arma::uvec::iterator ind_id = row_order.begin(); ind_id != ind_end;
       ++ind_id) {
    int indiv_idx = *ind_id;
    for (int iter_metro = 0; iter_metro < niter_metropolis; ++iter_metro) {
      // Sample the variance of the normal
      double rho_iter =
          rho * Rcpp::RcppArmadillo::sample(rho_values, 1, false)(0);

      arma::rowvec Pi_candidate = running_P.row(indiv_idx);
      arma::rowvec noise =
          arma::randn<arma::rowvec>(P.n_cols) * std::sqrt(rho_iter);

      Pi_candidate += noise;

      arma::rowvec Pi_old = running_P.row(indiv_idx);

      double log_u = log(Rcpp::runif(1)(0));

      arma::uvec Zi_vec = find(Z.row(indiv_idx) == 1, 1, "first");
      arma::uword Zi = Zi_vec(0);

      if (DEBUG_SAMPLE) {
        Rcpp::Rcout << "Z_" << indiv_idx + 1 << " = " << Zi + 1 << std::endl;
      }

      arma::rowvec ind_mu =
          mean_of_Pi_given_P_min_i_sigma(running_P, Sigma, sigma2, indiv_idx);
      arma::mat ind_cov =
          cov_of_Pi_given_P_min_i_sigma(running_P, Sigma, sigma2, indiv_idx);

      double log_accept = ilrInv_cpp(Pi_candidate, basis, true)(Zi) -
                          ilrInv_cpp(Pi_old, basis, true)(Zi);

      if (DEBUG_SAMPLE) {
        Rcpp::Rcout << "(multinom) accept_prob = " << exp(log_accept)
                    << std::endl;
      }

      log_accept += dmvnorm(Pi_candidate, ind_mu, ind_cov, true)(0) -
                    dmvnorm(Pi_old, ind_mu, ind_cov, true)(0);
      if (DEBUG_SAMPLE) {
        Rcpp::Rcout << "(multinom+mvnorm) accept_prob = " << exp(log_accept)
                    << " | u = " << exp(log_u) << std::endl;
      }

      if (log_u < log_accept) {
        if (DEBUG_SAMPLE) {
          Rcpp::Rcout << "Accepted !" << std::endl << "Old Pi = ";
          Pi_old.print(Rcpp::Rcout);
          Rcpp::Rcout << "New Pi = ";
          Pi_candidate.print(Rcpp::Rcout);
        }

        running_P.row(indiv_idx) = Pi_candidate;
        accepted_count++;
      }
    }
  }

  if (DEBUG_SAMPLE) {
    Rcpp::Rcout << "Mean accepted rate = "
                << (float)accepted_count / (float)(niter_metropolis * n)
                << std::endl;
  }

  return running_P;
};

//' A Metropolis-Hastings sampler for P with the clever proposition (C++)
//'
//' Proposes each row \eqn{P_i} from its full conditional
//' \eqn{N(\mu_i, \Sigma_i)} (as the R \code{sample_P_metropolis_trick}):
//' the Gaussian terms cancel in the acceptance ratio, leaving only the
//' multinomial prior ratio. This is the C++ version of
//' \code{sample_P_metropolis_trick}.
//' @param rho unused double for compatibility only
//' @inheritParams sample_P_metropolis_trick
//' @return an updated matrix of latent positions, of the same size
//' as \code{P}
//' @seealso [sample_P_metropolis_classical()],
//' [sample_P_metropolis_classical_cpp()],
//' [sample_P_metropolis_trick()]
//' @export
// [[Rcpp::export]]
arma::mat sample_P_metropolis_trick_cpp(arma::mat &P, arma::mat &Z,
                                        arma::mat &Sigma, double sigma2,
                                        bool minibatch = true,
                                        int niter_metropolis = 50,
                                        double rho = 1.0) {
  Rcpp::RNGScope scope;
  (void)rho;
  mat running_P = P;
  int n = running_P.n_rows;
  int K = running_P.n_cols + 1;
  arma::mat basis = default_Psi_function_cpp(K);

  arma::uvec row_order;

  if (minibatch) {
    row_order = arma::randperm(n);
  } else {
    row_order = arma::regspace<uvec>(0, n - 1);
  }

  // Loop over the individuals
  int accepted_count = 0;
  arma::uvec::iterator ind_end = row_order.end();
  for (arma::uvec::iterator ind_id = row_order.begin(); ind_id != ind_end;
       ++ind_id) {
    int indiv_idx = *ind_id;
    for (int iter_metro = 0; iter_metro < niter_metropolis; ++iter_metro) {
      // Conditional distribution of Pi given P_minus_i
      arma::rowvec ind_mu =
          mean_of_Pi_given_P_min_i_sigma(running_P, Sigma, sigma2, indiv_idx);
      arma::mat ind_cov =
          cov_of_Pi_given_P_min_i_sigma(running_P, Sigma, sigma2, indiv_idx);

      // Proposal drawn from the conditional Gaussian (sample_Pi_given): the
      // mvnorm terms cancel in the acceptance ratio, leaving only the
      // multinomial prior ratio.
      arma::rowvec Pi_old = running_P.row(indiv_idx);
      arma::rowvec Pi_candidate = rmvnorm(1, ind_mu, ind_cov).row(0);

      double log_u = log(Rcpp::runif(1)(0));

      arma::uvec Zi_vec = find(Z.row(indiv_idx) == 1, 1, "first");
      arma::uword Zi = Zi_vec(0);

      double log_accept = ilrInv_cpp(Pi_candidate, basis, true)(Zi) -
                          ilrInv_cpp(Pi_old, basis, true)(Zi);

      if (DEBUG_SAMPLE) {
        Rcpp::Rcout << "(multinom) accept_prob = " << exp(log_accept)
                    << " | u = " << exp(log_u) << std::endl;
      }

      if (log_u < log_accept) {
        if (DEBUG_SAMPLE) {
          Rcpp::Rcout << "Accepted !" << std::endl << "Old Pi = ";
          Pi_old.print(Rcpp::Rcout);
          Rcpp::Rcout << "New Pi = ";
          Pi_candidate.print(Rcpp::Rcout);
        }

        running_P.row(indiv_idx) = Pi_candidate;
        accepted_count++;
      }
    }
  }

  if (DEBUG_SAMPLE) {
    Rcpp::Rcout << "Mean accepted rate = "
                << (float)accepted_count / (float)(niter_metropolis * n)
                << std::endl;
  }

  return running_P;
};

//' A Metropolis-Hastings sampler for P with the clever proposition and block
//' update (C++)
//' @inheritParams sample_P_metropolis_trick_cpp
//' @param block_size the size of blocks to update simultaneously (default: 1
//' for individual updates) ' @return an updated matrix of latent positions, of
//' the same size ' as \code{P} ' @seealso [sample_P_metropolis_classical()], '
//' [sample_P_metropolis_classical_cpp()], ' [sample_P_metropolis_trick()] '
//' @export
// [[Rcpp::export]]
arma::mat sample_P_metropolis_trick_cpp_block(arma::mat &P, arma::mat &Z,
                                              arma::mat &Sigma, double sigma2,
                                              bool minibatch = true,
                                              int niter_metropolis = 50,
                                              double rho = 1.0,
                                              int block_size = 1) {
  Rcpp::RNGScope scope;
  (void)rho;
  mat running_P = P;
  int n = running_P.n_rows;
  int K = running_P.n_cols + 1;
  arma::mat basis = default_Psi_function_cpp(K);

  arma::uvec row_order;

  if (minibatch) {
    row_order = arma::randperm(n);
  } else {
    row_order = arma::regspace<uvec>(0, n - 1);
  }

  // Loop over the individuals in blocks
  int accepted_count = 0;
  int total_proposals = 0;

  // Process individuals in blocks
  arma::uvec::iterator ind_end = row_order.end();
  arma::uvec::iterator ind_id = row_order.begin();

  while (ind_id != ind_end) {
    // Determine the current block size (may be smaller for the last block)
    int current_block_size = std::min(block_size, (int)(ind_end - ind_id));

    // Collect indices for the current block
    arma::uvec block_indices(current_block_size);
    for (int i = 0; i < current_block_size; ++i) {
      block_indices(i) = *(ind_id + i);
    }

    // For each iteration of Metropolis-Hastings within the block
    for (int iter_metro = 0; iter_metro < niter_metropolis; ++iter_metro) {
      // Store old values for all individuals in the block
      arma::mat old_P_block = running_P.rows(block_indices);

      // Generate proposals for all individuals in the block
      arma::mat new_P_block(old_P_block.n_rows, old_P_block.n_cols);

      // For each individual in the block, compute conditional distribution and
      // propose
      for (int i = 0; i < current_block_size; ++i) {
        int indiv_idx = block_indices(i);

        // Conditional distribution of Pi given P_minus_i
        arma::rowvec ind_mu =
            mean_of_Pi_given_P_min_i_sigma(running_P, Sigma, sigma2, indiv_idx);
        arma::mat ind_cov =
            cov_of_Pi_given_P_min_i_sigma(running_P, Sigma, sigma2, indiv_idx);

        // Proposal drawn from the conditional Gaussian
        arma::rowvec Pi_old = running_P.row(indiv_idx);
        arma::rowvec Pi_candidate = rmvnorm(1, ind_mu, ind_cov).row(0);
        new_P_block.row(i) = Pi_candidate;
      }

      // Compute acceptance probability for the entire block
      double log_accept_total = 0.0;
      for (int i = 0; i < current_block_size; ++i) {
        int indiv_idx = block_indices(i);
        arma::rowvec Pi_old = old_P_block.row(i);
        arma::rowvec Pi_candidate = new_P_block.row(i);

        arma::uvec Zi_vec = find(Z.row(indiv_idx) == 1, 1, "first");
        arma::uword Zi = Zi_vec(0);

        double log_accept_single = ilrInv_cpp(Pi_candidate, basis, true)(Zi) -
                                   ilrInv_cpp(Pi_old, basis, true)(Zi);

        log_accept_total += log_accept_single;
      }

      double log_u = log(Rcpp::runif(1)(0));
      total_proposals += current_block_size;

      if (log_u < log_accept_total) {
        // Accept the entire block
        running_P.rows(block_indices) = new_P_block;
        accepted_count += current_block_size;
      }
    }

    // Move to the next block
    std::advance(ind_id, current_block_size);
  }

  if (DEBUG_SAMPLE) {
    Rcpp::Rcout << "Mean accepted rate = "
                << (float)accepted_count / (float)total_proposals << std::endl;
  }

  return running_P;
};

//' A fast Metropolis sampler for P with the clever proposition (C++)
//'
//' Same proposal as [sample_P_metropolis_trick_cpp()] but computed from the
//' precision (inverse covariance) matrix `Theta` instead of `Sigma`. For each
//' row \eqn{i} the full conditional is \eqn{N(\mu_i, \Sigma_i)} with
//' \deqn{\mu_i = -\frac{1}{\Theta_{ii}} \Theta_{i,-i} P_{-i}, \qquad
//' \Sigma_i = \frac{\sigma^2}{\Theta_{ii}} I,}
//' so the conditional mean is obtained with a single matrix-vector product
//' (\eqn{O(nK)}) instead of a fresh \eqn{O(n^3)} solve of `Sigma[-i,-i]` per
//' row. As in the ``trick'' samplers, drawing the proposal from this full
//' conditional makes the Gaussian terms cancel in the acceptance ratio,
//' leaving only the multinomial prior ratio.
//'
//' @param P a matrix of size \eqn{n_1 \times K-1} specifying latent positions
//' @param Z a matrix of size \eqn{n_1\times K} with a single 1 per line
//'   indicating the membership of row node \eqn{i}
//' @param Theta the precision matrix, i.e. \eqn{\Sigma^{-1}}
//' @param sigma2 a variance parameter indicating the variance between the
//'   K-1 columns of P
//' @param minibatch a boolean indicating wether to update the rows in the
//'   lexical order or to sample at each iteration. Default to TRUE
//' @param niter_metropolis an integer specifying the number of metropolis
//'   iterations to perform. Defaults to 50.
//' @return an updated matrix of latent positions, of the same size as
//'   \code{P}
// [[Rcpp::export(name = "sample_P_metropolis_trick_cpp_fast_impl")]]
arma::mat sample_P_metropolis_trick_cpp_fast(arma::mat &P, arma::mat &Z,
                                             arma::mat &Theta, double sigma2,
                                             bool minibatch = true,
                                             int niter_metropolis = 50) {
  Rcpp::RNGScope scope;
  mat running_P = P;
  int n = running_P.n_rows;
  int q = running_P.n_cols;
  int K = q + 1;
  mat basis = default_Psi_function_cpp(K);

  uvec row_order;
  if (minibatch) {
    row_order = random_perm_r(n);
  } else {
    row_order = regspace<uvec>(0, n - 1);
  }

  int accepted_count = 0;
  for (uvec::iterator ind_id = row_order.begin(); ind_id != row_order.end();
       ++ind_id) {
    int i = *ind_id;
    double Theta_ii = Theta(i, i);
    if (Theta_ii <= 0) {
      Rcpp::stop("Theta must be positive definite: diagonal entry %d <= 0.",
                 i + 1);
    }
    double sd = std::sqrt(sigma2 / Theta_ii);
    // Precision-matrix row products (Theta only heterogeneous through row i).
    arma::subview_row<double> Theta_i = Theta.row(i);
    // Conditional mean: mu = -(1/Theta_ii) * Theta[i,-i] %*% P[-i,]
    rowvec mu(q, arma::fill::zeros);
    for (int j = 0; j < n; ++j) {
      if (j == i) continue;
      mu -= (Theta_i(j) / Theta_ii) * running_P.row(j);
    }

    for (int iter_metro = 0; iter_metro < niter_metropolis; ++iter_metro) {
      // Draw along R's RNG stream exactly as rmvnorm does in the Sigma-based
      // trick samplers (one Rcpp::rnorm per coordinate), so that the two share
      // the same RNG consumption and the same `.Random.seed` evolution.
      rowvec Pi_candidate(q);
      for (int k = 0; k < q; ++k) {
        Pi_candidate(k) = Rcpp::rnorm(1, 0.0, 1.0)(0) * sd + mu(k);
      }
      rowvec Pi_old = running_P.row(i);

      double log_u = log(Rcpp::runif(1)(0));

      uvec Zi_vec = find(Z.row(i) == 1, 1, "first");
      uword Zi = Zi_vec(0);

      double log_accept = ilrInv_cpp(Pi_candidate, basis, true)(Zi) -
                          ilrInv_cpp(Pi_old, basis, true)(Zi);

      if (DEBUG_SAMPLE) {
        Rcpp::Rcout << "(multinom) accept_prob = " << exp(log_accept)
                    << " | u = " << exp(log_u) << std::endl;
      }

      if (log_u < log_accept) {
        if (DEBUG_SAMPLE) {
          Rcpp::Rcout << "Accepted !" << std::endl << "Old Pi = ";
          Pi_old.print(Rcpp::Rcout);
          Rcpp::Rcout << "New Pi = ";
          Pi_candidate.print(Rcpp::Rcout);
        }
        running_P.row(i) = Pi_candidate;
        accepted_count++;
      }
    }
  }

  if (DEBUG_SAMPLE) {
    Rcpp::Rcout << "Mean accepted rate = "
                << (float)accepted_count / (float)(niter_metropolis * n)
                << std::endl;
  }

  return running_P;
};

//' A fast Metropolis sampler for P under soft memberships (C++)
//'
//' Proposes every row of `P` from its full conditional given the precision
//' matrix `Theta` (see [sample_P_metropolis_trick_cpp_fast()]) and accepts
//' with the *mixed-membership* ratio
//' \deqn{\sum_k Z^\text{soft}_{i,k}\,
//' \left[\log p_k(\Pi_i^\text{cand}) - \log p_k(\Pi_i)\right],}
//' where \eqn{p_k} are the probabilities obtained through the inverse pivot
//' coordinates. This is the fast, precision-based analogue of
//' [sample_P_metropolis_trick_soft()].
//'
//' @param P a matrix of size \eqn{n_1 \times K-1} specifying latent positions
//' @param Z a matrix (\eqn{n_1 \times K}) of **soft** membership
//'   probabilities (rows sum to 1)
//' @param Theta the precision matrix, i.e. \eqn{\Sigma^{-1}}
//' @param sigma2 a variance parameter indicating the variance between the
//'   K-1 columns of P
//' @param minibatch a boolean indicating wether to update the rows in the
//'   lexical order or to sample at each iteration. Default to TRUE
//' @param niter_metropolis an integer specifying the number of metropolis
//'   iterations to perform. Defaults to 50.
//' @return an updated matrix of latent positions, of the same size as
//'   \code{P}
// [[Rcpp::export(name = "sample_P_metropolis_trick_soft_cpp_fast_impl")]]
arma::mat sample_P_metropolis_trick_soft_cpp_fast(arma::mat &P, arma::mat &Z,
                                                  arma::mat &Theta,
                                                  double sigma2,
                                                  bool minibatch = true,
                                                  int niter_metropolis = 50) {
  Rcpp::RNGScope scope;
  mat running_P = P;
  int n = running_P.n_rows;
  int q = running_P.n_cols;
  int K = q + 1;
  mat basis = default_Psi_function_cpp(K);

  uvec row_order;
  if (minibatch) {
    row_order = random_perm_r(n);
  } else {
    row_order = regspace<uvec>(0, n - 1);
  }

  int accepted_count = 0;
  for (uvec::iterator ind_id = row_order.begin(); ind_id != row_order.end();
       ++ind_id) {
    int i = *ind_id;
    double Theta_ii = Theta(i, i);
    if (Theta_ii <= 0) {
      Rcpp::stop("Theta must be positive definite: diagonal entry %d <= 0.",
                 i + 1);
    }
    double sd = std::sqrt(sigma2 / Theta_ii);
    arma::subview_row<double> Theta_i = Theta.row(i);
    rowvec mu(q, arma::fill::zeros);
    for (int j = 0; j < n; ++j) {
      if (j == i) continue;
      mu -= (Theta_i(j) / Theta_ii) * running_P.row(j);
    }

    for (int iter_metro = 0; iter_metro < niter_metropolis; ++iter_metro) {
      // Draw along R's RNG stream exactly as rmvnorm does in the Sigma-based
      // trick samplers (one Rcpp::rnorm per coordinate), so that the two share
      // the same RNG consumption and the same `.Random.seed` evolution.
      rowvec Pi_candidate(q);
      for (int k = 0; k < q; ++k) {
        Pi_candidate(k) = Rcpp::rnorm(1, 0.0, 1.0)(0) * sd + mu(k);
      }
      rowvec Pi_old = running_P.row(i);

      rowvec log_pi_cand = ilrInv_cpp(Pi_candidate, basis, true);
      rowvec log_pi_old = ilrInv_cpp(Pi_old, basis, true);
      double log_accept = dot(Z.row(i), log_pi_cand - log_pi_old);

      double log_u = log(Rcpp::runif(1)(0));

      if (DEBUG_SAMPLE) {
        Rcpp::Rcout << "(soft multinom) accept_prob = " << exp(log_accept)
                    << " | u = " << exp(log_u) << std::endl;
      }

      if (log_u < log_accept) {
        if (DEBUG_SAMPLE) {
          Rcpp::Rcout << "Accepted !" << std::endl << "Old Pi = ";
          Pi_old.print(Rcpp::Rcout);
          Rcpp::Rcout << "New Pi = ";
          Pi_candidate.print(Rcpp::Rcout);
        }
        running_P.row(i) = Pi_candidate;
        accepted_count++;
      }
    }
  }

  if (DEBUG_SAMPLE) {
    Rcpp::Rcout << "Mean accepted rate = "
                << (float)accepted_count / (float)(niter_metropolis * n)
                << std::endl;
  }

  return running_P;
};
