// Regularized Student-t
data {
  int<lower = 0> n_stan; // num obs
  int<lower = 0> p_stan; // number of covariates
  real y_stan[n_stan]; // outcome
  matrix[n_stan, p_stan] x_standardized_stan; //covariates (no intercept)
  real<lower = 0> dof_stan;
  real<lower = 0> beta_scale_stan;
  // d in Section 2 is InvGamma(slab_dof_stan/2,slab_scale_stan^2*slab_dof_stan/2)
  real<lower = 0> slab_dof_stan; //
  real<lower = 0> slab_scale_stan;
  real<lower = 0> mu_sd_stan; // prior standard deviation on main intercept
  int<lower = 0,upper = 1> only_prior;//if 1, ignore the model and data
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
  vector[p_stan] beta_raw;
  vector<lower = 0>[p_stan] lambda_scale_sq;//shrinkage factor for betas
  real<lower = 0> slab; // d^2, Section 2
  real<lower = 0> sigma;
}
transformed parameters {
  vector[p_stan] beta;//shrunk regression coefficients
  vector<lower = 0>[p_stan] theta;
  real<lower = 0> slab_copy;// to gracefully allow for infinite slab_dof_stan
  if(!is_inf(slab_dof_stan)) {
    slab_copy = slab;
  } else {
    slab_copy = slab_scale_stan^2;
  }
  theta = 1 ./ sqrt(1 / slab_copy + (1 ./ (beta_scale_stan^2 * sigma^2 * lambda_scale_sq)));
  beta = theta .* beta_raw;
}
model {
  beta_raw ~ normal(0.0, 1.0);
  lambda_scale_sq ~ inv_gamma(dof_stan/2.0, dof_stan/2.0);
  slab ~ inv_gamma(slab_ig_shape, slab_ig_scale);
  sigma ~ student_t(1, 0.0, 5.0);
  mu ~ logistic(0.0, mu_sd_stan);
  if(only_prior == 0)
  y_stan ~ normal(mu + x_standardized_stan * beta, sigma);
}
