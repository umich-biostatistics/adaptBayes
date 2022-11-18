// Regularized Student-t
data {
  int<lower = 0> n_stan; // num obs
  int<lower = 0> p_stan; // number of covariates
  int<lower = 0,upper=1> y_stan[n_stan]; // outcome
  matrix[n_stan, p_stan] x_standardized_stan; //covariates (no intercept)
  real<lower = 0> dof_stan;
  real<lower = 0> beta_scale_stan;
  real<lower = 0> slab_precision_stan; // 1/sqrt(slab_precision_stan) is the maximum possible scale
  real<lower = 0> mu_sd_stan; // prior standard deviation on main intercept
  int<lower = 0,upper = 1> only_prior;//if 1, ignore the model and data
}
parameters {
  real mu;
  vector[p_stan] beta_raw;
  vector<lower = 0>[p_stan] lambda_scale_sq;//shrinkage factor for betas
  }
transformed parameters {
  vector[p_stan] beta;//shrunk regression coefficients
  vector<lower = 0,upper = sqrt(1/slab_precision_stan)>[p_stan] theta;
  theta = 1 ./ sqrt(slab_precision_stan + (1 ./ (beta_scale_stan^2 * lambda_scale_sq)));
  beta = theta .* beta_raw;
}
model {
  beta_raw ~ normal(0.0, 1.0);
  // We multiply the scale parameter (the second parameter) to take into account
  // the contribution from "sigma", which
  // isn't really a parameter in a logistic glm but follows the maximum variance
  // assumption that Pirronen et al. suggest in their discussion on the choice
  // of the scale parameter
  lambda_scale_sq ~ inv_gamma(dof_stan/2.0, dof_stan);
  mu ~ logistic(0.0, mu_sd_stan);
  if(only_prior == 0)
    y_stan ~ bernoulli_logit(mu + x_standardized_stan * beta);
}
