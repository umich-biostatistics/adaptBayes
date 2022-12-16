//Naive Adaptive Bayes (stable version)
// In this simple version psi and phi are identically equal to 1.
// Refer to nab_binomial.stan for complete version
data {
  int<lower = 1> n_stan; // n_curr
  int<lower = 1> p_stan; // number of original covariates
  int<lower = 0> q_stan; // number of augmented (or added) covariates
  int<lower = 0, upper = 1> y_stan[n_stan]; // outcome
  matrix[n_stan, p_stan + q_stan] x_standardized_stan; // covariates (no intercept)
  vector[p_stan] alpha_prior_mean_stan; // m_alpha
  matrix[p_stan, p_stan] alpha_prior_cov_stan; // S_alpha
  vector[p_stan] sqrt_eigenval_hist_var_stan; // sqrt of eigenvalues of S_alpha; D^1/2 in Equation (S5)
  matrix[p_stan, p_stan] eigenvec_hist_var_stan; // eigenvectors of S_alpha; Q^T in Equation (S5)
  real<lower = 0> local_dof_stan; // dof of pi(lambda), = 1
  real<lower = 0> global_dof_stan; // dof of pi(tau), = 1
  real<lower = 0> beta_orig_scale_stan; // c, Section 2
  real<lower = 0> beta_aug_scale_stan; // c, Section 2
  real<lower = 0> beta_aug_scale_tilde_stan; // tilde~c = 0.05, Section 4.1
  // d in Section 2 is InvGamma(slab_dof_stan/2,slab_scale_stan^2*slab_dof_stan/2)
  real<lower = 0> slab_dof_stan; //
  real<lower = 0> slab_scale_stan;
  vector<lower = 0>[p_stan] scale_to_variance225; // not used in simple version but needed for compatibility with glm_sab
  int<lower = 0, upper = 1> phi_prior_type; // not used in simple version but needed for compatibility with glm_sab
  real<lower = 0, upper = 1> phi_mean_stan; // not used in simple version but needed for compatibility with glm_sab
  real<lower = 0> phi_sd_stan; // not used in simple version but needed for compatibility with glm_sab
  real<lower = 0> eta_param_stan; // not used in simple version but needed for compatibility with glm_sab
  real<lower = 0> mu_sd_stan; // prior standard deviation on intercept
  int<lower = 0, upper = 1> only_prior;//if 1, ignore the model and data and generate from the prior only
}
transformed data {
  vector[p_stan] zero_vec;
  real<lower = 0> slab_ig_shape;
  real<lower = 0> slab_ig_scale;
  zero_vec = rep_vector(0.0, p_stan);
  if(!is_inf(slab_dof_stan)) {
    slab_ig_shape = slab_dof_stan / 2.0;
    slab_ig_scale = slab_scale_stan^2 * slab_dof_stan / 2.0;
  } else {
    // slab is not used in the model here but we use a proper prior to avoid
    // numerical issues
    slab_ig_shape = 2.5;
    slab_ig_scale = 2.5;
  }
}
parameters {
  real mu;
  vector[p_stan] beta_raw_orig; // unscaled
  vector[q_stan] beta_raw_aug; // unscaled
  vector<lower = 0>[q_stan] lambda_aug_tilde_sq; // tilde lambda^2
  real<lower = 0> tau_glob;// tau
  vector<lower = 0>[p_stan] lambda_orig;// lambda^o
  vector<lower = 0>[q_stan] lambda_aug;// lambda^a
  real<lower = 0> slab; // d^2, Section 2
}
transformed parameters {
  vector[p_stan] normalized_beta; //entire LHS of Equation (S5)
  vector[p_stan] beta_orig; // scaled
  vector[q_stan] beta_aug; // scaled
  vector[p_stan + q_stan] beta;
  vector<lower = 0>[p_stan] theta_orig; // theta
  vector<lower = 0>[q_stan] theta_aug; // theta
  matrix[p_stan,p_stan] normalizing_cov;// S_alpha  + Theta^o
  real<lower = 0> slab_copy;// to gracefully allow for infinite slab_dof_stan
  if(!is_inf(slab_dof_stan)) {
    slab_copy = slab;
  } else {
    slab_copy = slab_scale_stan^2;
  }
  theta_orig = 1 ./ sqrt(1 / slab_copy + zero_vec);
  theta_aug = 1 ./ sqrt(1 / slab_copy + (1.0 ./ (beta_aug_scale_tilde_stan^2 * lambda_aug_tilde_sq)));
  beta_orig = theta_orig .* beta_raw_orig;
  beta_aug = theta_aug .* beta_raw_aug;
  beta = append_row(beta_orig, beta_aug);
  normalizing_cov = alpha_prior_cov_stan;
  for(i in 1:p_stan) {
    normalizing_cov[i,i] = normalizing_cov[i,i] + theta_orig[i]^2;
  }
  normalized_beta = eigenvec_hist_var_stan * (beta_orig - alpha_prior_mean_stan);
}
model {
  beta_raw_orig ~ normal(0.0, 1.0);
  beta_raw_aug ~ normal(0.0, 1.0);
  lambda_aug_tilde_sq ~ inv_gamma(0.5, 0.5);
  tau_glob ~ student_t(global_dof_stan, 0.0, 1.0);
  // The 2.0 scaler below represents the contribution from "sigma", which
  // isn't really a parameter in a logistic glm but follows the maximum variance
  // assumption that Pirronen et al. suggest in their discussion on the choice
  // of the scale parameter
  lambda_orig ~ student_t(local_dof_stan, 0.0, 2.0 * beta_orig_scale_stan);
  lambda_aug ~ student_t(local_dof_stan, 0.0, 2.0 * beta_aug_scale_stan);
  slab ~ inv_gamma(slab_ig_shape, slab_ig_scale);
  mu ~ logistic(0.0, mu_sd_stan);
  // Equation (S5) This is the naive adaptive prior contribution
  normalized_beta ~ normal(0.0, sqrt_eigenval_hist_var_stan);
  // Z_NAB (Normalizing constant)
  target += -(1.0 * multi_normal_lpdf(alpha_prior_mean_stan|zero_vec, normalizing_cov));
  if(only_prior == 0) {
    y_stan ~ bernoulli_logit_glm(x_standardized_stan, mu, beta);
  }
}
