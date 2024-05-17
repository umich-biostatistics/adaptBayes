// Regularized Horseshoe Prior (stable version)
data {
  int<lower = 1> n_stan; // n_curr
  int<lower = 1> p_stan; // number of original covariates
  int<lower = 0> q_stan; // number of augmented (or added) covariates
  array[n_stan] real y_stan; // outcome
  matrix[n_stan, p_stan + q_stan] x_standardized_stan; //covariates (no intercept)
  real<lower = 0> local_dof_stan; // dof of pi(lambda), = 1
  real<lower = 0> global_dof_stan; // dof of pi(tau), = 1
  real<lower = 0> beta_orig_scale_stan; // c, Section 2
  real<lower = 0> beta_aug_scale_stan; // c, Section 2
  // d in Section 2 is InvGamma(slab_dof_stan/2,slab_scale_stan^2*slab_dof_stan/2)
  real<lower = 0> slab_dof_stan; //
  real<lower = 0> slab_scale_stan;
  real<lower = 0> mu_sd_stan; // prior standard deviation on main intercept
  vector<lower = 0,upper = 1>[n_stan] intercept_offset_stan;//allow for two intercepts
  int<lower = 0,upper = 1> only_prior;//if 1, ignore the model and data and generate from the prior only
}
transformed data {
  real<lower = 0> slab_ig_shape;
  real<lower = 0> slab_ig_scale;
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
  real mu_offset;
  vector[p_stan + q_stan] beta_raw;
  real<lower = 0> tau_glob; // tau
  vector<lower = 0>[p_stan] lambda_orig;// lambda^o
  vector<lower = 0>[q_stan] lambda_aug; // lambda^a
  real<lower = 0> slab; // d^2, Section 2
  real<lower = 0> sigma;
}
transformed parameters {
  vector[p_stan + q_stan] beta;
  vector<lower = 0>[p_stan] theta_orig;//theta
  vector<lower = 0>[q_stan] theta_aug;//theta
  real<lower = 0> slab_copy;// to gracefully allow for infinite slab_dof_stan
  if(!is_inf(slab_dof_stan)) {
    slab_copy = slab;
  } else {
    slab_copy = slab_scale_stan^2;
  }
  theta_orig = 1 ./ sqrt(1 / slab_copy + (1 ./ (tau_glob^2 * sigma^2 * square(lambda_orig))));
  theta_aug = 1 ./ sqrt(1 / slab_copy + (1 ./ (tau_glob^2 * sigma^2 * square(lambda_aug))));
  beta = append_row(theta_orig, theta_aug) .* beta_raw;
}
model {
  beta_raw ~ normal(0.0, 1.0);
  tau_glob ~ student_t(global_dof_stan, 0.0, 1.0);
  lambda_orig ~ student_t(local_dof_stan, 0.0, beta_orig_scale_stan);
  lambda_aug ~ student_t(local_dof_stan, 0.0, beta_aug_scale_stan);
  slab ~ inv_gamma(slab_ig_shape, slab_ig_scale);
  sigma ~ student_t(1, 0.0, 5.0);
  mu ~ logistic(0.0, mu_sd_stan);
  mu_offset ~ logistic(0.0, mu_sd_stan / 2);
  if(only_prior == 0) {
    y_stan ~ normal(mu + intercept_offset_stan * mu_offset + x_standardized_stan * beta, sigma);
  }
}
