//Sensible Adaptive Bayes (stable version)
// This is only (Equation 3.9) from the manuscript, i.e. there is no
// horseshoe shrinkage
data {
  int<lower = 1> n_stan; // n_curr
  int<lower = 1> p_stan; // number of original covariates
  int<lower = 0> q_stan; // number of augmented (or added) covariates
  real y_stan[n_stan]; // outcome
  matrix[n_stan, p_stan + q_stan] x_standardized_stan; // covariates (no intercept)
  matrix[p_stan, q_stan] aug_projection_stan; // {v^o}^(-1) %*% (E[V^a|v^o]- E[V^a|v^o = 0]), Eqn (9)
  vector[p_stan] alpha_prior_mean_stan; // m_alpha
  matrix[p_stan, p_stan] alpha_prior_cov_stan; // S_alpha
  vector[p_stan] sqrt_eigenval_hist_var_stan; // sqrt of eigenvalues of S_alpha; D^1/2 in Equation (S6)
  matrix[p_stan, p_stan] eigenvec_hist_var_stan; // eigenvectors of S_alpha; Q^T in Equation (S6)
  vector<lower = 0>[p_stan] scale_to_variance225; //Equation (S6); equal to diag(S_alpha) / 225;
  int<lower = 0, upper = 1> phi_prior_type; //if 1, phi is apriori truncated normal; if 0, phi is apriori beta
  real<lower = 0, upper = 1> phi_mean_stan; //
  real<lower = 0> phi_sd_stan; //
  real<lower = 0> eta_param_stan; // eta ~ IG(eta_param_stan, eta_param_stan)
  real<lower = 0> mu_sd_stan; // prior standard deviation on intercept
  int<lower = 0, upper = 1> only_prior;//if 1, ignore the model and data and generate from the prior only
}
transformed data {
  vector[p_stan] zero_vec;
  real phi_beta_shape1;
  real phi_beta_shape2;
  real phi_mean_stan_trunc;
  real phi_sd_stan_trunc;
  zero_vec = rep_vector(0.0, p_stan);
  // These next four lines will only be used if phi is beta distributed
  // Depending on the values of phi_mean_stan or phi_sd_stan, the resulting
  // shape parameters could give stan fits
  phi_mean_stan_trunc = fmax(1e-6, fmin(1 - 1e-6, phi_mean_stan));
  phi_sd_stan_trunc = fmax(1e-6, fmin(sqrt(phi_mean_stan_trunc * (1 - phi_mean_stan_trunc)) - 1e-6, phi_sd_stan));
  phi_beta_shape1 = phi_mean_stan_trunc * (phi_mean_stan_trunc * (1 - phi_mean_stan_trunc) / phi_sd_stan_trunc^2 - 1);
  phi_beta_shape2 = (1 - phi_mean_stan_trunc) * (phi_mean_stan_trunc * (1 - phi_mean_stan_trunc) / phi_sd_stan_trunc^2 - 1);
}
parameters {
  real mu;
  vector[p_stan] beta_orig; //
  vector[q_stan] beta_aug; //
  real<lower = 0> eta;
  real<lower = 0, upper = 1> phi;// phi
  real<lower = 0> sigma;
}
transformed parameters {
  vector[p_stan] normalized_beta; //entire LHS of Equation (S6)
  vector[p_stan + q_stan] beta;
  vector<lower = 0>[p_stan] hist_orig_scale;//Diagonal of Gamma^{-1} in LHS of Eqn (S6)
  real<lower = 0, upper = 1> phi_copy;// to gracefully allow for zero-valued prior standard deviation
  real<lower = 0> eta_copy;// to gracefully allow for point mass on 1
  if(phi_sd_stan > 0) {
    phi_copy = phi;
  } else {
    phi_copy = phi_mean_stan;
  }
  if(!is_inf(eta_param_stan)) {
    eta_copy = eta;
  } else {
    eta_copy = 1;
  }
  beta = append_row(beta_orig, beta_aug);
  hist_orig_scale = 1 ./ sqrt(scale_to_variance225 * (1 - phi_copy) + phi_copy / eta_copy);
  normalized_beta = eigenvec_hist_var_stan * (beta_orig + aug_projection_stan * beta_aug - alpha_prior_mean_stan) ./ hist_orig_scale;
}
model {
  if(!is_inf(eta_param_stan)) {
    eta ~ inv_gamma(eta_param_stan, eta_param_stan);
  } else {
    // eta is not used in this case but we need a proper prior to avoid sampling issues
    eta ~ inv_gamma(1.0, 1.0);
  }
  mu ~ logistic(0.0, mu_sd_stan);
  if(phi_sd_stan > 0 && phi_prior_type == 1) {
    phi ~ normal(phi_mean_stan, phi_sd_stan);
  } else if(phi_sd_stan > 0 && phi_prior_type == 0) {
    phi ~ beta(phi_beta_shape1, phi_beta_shape2);
  }
  sigma ~ student_t(1, 0.0, 5.0);
  // Equation (S6) The next two lines together comprise the sensible adaptive prior contribution
  normalized_beta ~ normal(0.0, sqrt_eigenval_hist_var_stan);
  // Scaling normalized_beta to be independent ends up dropping a necessary determinant calculation: we add that back in here:
  target += -(1.0 * sum(log(hist_orig_scale)));
  if(only_prior == 0) {
    y_stan ~ normal(mu + x_standardized_stan * beta, sigma);
  }
}
