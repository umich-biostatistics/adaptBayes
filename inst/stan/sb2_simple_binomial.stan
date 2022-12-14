//Sensible Bayes for binomial outcomes, version 2
// Version 2 uses the psi formulation
// This is only (Equation 3.9) from the manuscript, i.e. there is no
// horseshoe shrinkage
// In this simple version eta and phi are identically equal to 1.
// Refer to sb2_binomial.stan for complete version
data {
  int<lower = 1> n_stan; // n_curr
  int<lower = 1> p_stan; // number of original covariates
  int<lower = 0> q_stan; // number of augmented (or added) covariates
  int<lower = 0, upper = 1> y_stan[n_stan]; // outcome
  matrix[n_stan, p_stan + q_stan] x_standardized_stan; // covariates (no intercept)
  matrix[p_stan, q_stan] aug_projection_stan; // {v^o}^(-1) %*% (E[V^a|v^o]- E[V^a|v^o = 0]), Eqn (9)
  vector[p_stan] alpha_prior_mean_stan; // m_alpha
  matrix[p_stan, p_stan] alpha_prior_cov_stan; // S_alpha
  vector[p_stan] sqrt_eigenval_hist_var_stan; // sqrt of eigenvalues of S_alpha; D^1/2 in Equation (S6)
  matrix[p_stan, p_stan] eigenvec_hist_var_stan; // eigenvectors of S_alpha; Q^T in Equation (S6)
  vector<lower = 0>[p_stan] scale_to_variance225; // not used in simple version but needed for compatibility with glm_sb2
  int<lower = 0, upper = 1> phi_prior_type; // not used in simple version but needed for compatibility with glm_sb2
  real<lower = 0, upper = 1> phi_mean_stan; // not used in simple version but needed for compatibility with glm_sb2
  real<lower = 0> phi_sd_stan; // not used in simple version but needed for compatibility with glm_sb2
  real psi_mean_stan; // not used in simple version but needed for compatibility with glm_sb2
  real<lower = 0> psi_sd_stan; // not used in simple version but needed for compatibility with glm_sb2
  real<lower = 0> mu_sd_stan; // prior standard deviation on intercept
  int<lower = 0, upper = 1> only_prior;//if 1, ignore the model and data and generate from the prior only
}
transformed data {
  vector[p_stan] zero_vec;
  zero_vec = rep_vector(0.0, p_stan);
}
parameters {
  real mu;
  vector[p_stan] beta_orig; //
  vector[q_stan] beta_aug; //
}
transformed parameters {
  vector[p_stan] normalized_beta; //entire LHS of Equation (S6)
  vector[p_stan + q_stan] beta;
  beta = append_row(beta_orig, beta_aug);
  normalized_beta = eigenvec_hist_var_stan * (beta_orig + aug_projection_stan * beta_aug - alpha_prior_mean_stan);
}
model {
  mu ~ logistic(0.0, mu_sd_stan);
  // Equation (S6) This is the sensible adaptive prior contribution
  normalized_beta ~ normal(0.0, sqrt_eigenval_hist_var_stan);
  // Z_SAB (Normalizing constant)
  if(only_prior == 0) {
    y_stan ~ bernoulli_logit(mu + x_standardized_stan * beta);
  }
}
