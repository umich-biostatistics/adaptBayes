
#DESCRIPTION: Simulator function for drawing binary outcomes (historical, current, and new [for validation]),
#the probabilities of which are logistic-linear functions of normal and/or bernoulli  distributed
#predictors.
#
#
#ARGUMENTS:
#
#n_hist (pos. integer) size of historical data; n_hist in the paper
#
#n_curr (pos. integer) size of current data; n_curr in the paper
#
#n_new (pos. integer) size of testing data (for prediction)
#
#true_mu_hist (real) true intercept for generating historical model. mu_hist in Boonstra and Barbaro
#
#true_mu_curr (real) true intercept for generating current / new model. mu in Boonstra and Barbaro
#
#true_betas_orig (vector) true regression coefficients corresponding to original covariates. Beta^o in Boonstra and Barbaro
#
#true_betas_aug (vector) true regression coefficients corresponding to augmented covariates. Beta^a in Boonstra and Barbaro
#
#covariate_args (list) the named arguments are simple ways to govern the distribution of the predictors.
#x_correlation is the normal correlation between any pair of predictors; x_orig_binom is the
#integer indices (any subset of 1, ..., length(true_betas_orig)) indicating which of the original
#covariates should be transformed to binary values based upon being less than or greater than zero;
#x_aug_binom is the analogous set of indices for the augmented predictors (any subset of 1, ..., length(true_betas_aug))

draw_data = function(n_hist = 150,
                     n_curr = 50,
                     n_new = 100,
                     true_mu_hist = 0,
                     true_mu_curr = 0,
                     true_betas_orig = 0,
                     true_betas_aug = 0,
                     covariate_args = list(x_correlation = 0,
                                           x_orig_binom = NULL,
                                           x_aug_binom = NULL)
) {

  p = length(true_betas_orig);
  q = length(true_betas_aug);
  stopifnot(length(true_mu_hist) == 1 && length(true_mu_curr) == 1);

  x_all = matrix(rnorm((n_hist + n_curr + n_new)*(p+q)),nrow=n_hist + n_curr + n_new)%*%chol(diag(1 - covariate_args$x_correlation,p+q) + covariate_args$x_correlation);#original covariates are N(0,1)
  x_all_orig = x_all[,1:p,drop = F];
  #Binary covariates will be -1 or 1 with equal prevalence (assuming the latent normal is mean zero),
  #which will result in a random variable with mean zero and variance 1
  if(length(covariate_args$x_orig_binom)) {
    x_all_orig[,covariate_args$x_orig_binom] = 2 * (x_all_orig[,covariate_args$x_orig_binom,drop = F] > 0) - 1;
  }
  x_all_aug = x_all[,(p + 1):(p + q),drop = F];
  if(length(covariate_args$x_aug_binom)) {
    x_all_aug[,covariate_args$x_aug_binom] = 2 * (x_all_aug[,covariate_args$x_aug_binom,drop = F] > 0) - 1;
  }

  #Linear predictors (two for each observation: contribution from original covariates and augmented covariates)
  lin_pred_x_orig = drop(x_all_orig%*%true_betas_orig);
  lin_pred_x_aug = drop(x_all_aug%*%true_betas_aug);
  #Note the difference in intercepts between historical data and current data
  risk_all = 1/(1+exp(-c(rep(true_mu_hist,n_hist),rep(true_mu_curr,n_curr+n_new)) - lin_pred_x_orig - lin_pred_x_aug));

  y_all =  rbinom(n_hist + n_curr + n_new, 1, risk_all);

  x_hist_orig = x_all_orig[1:n_hist,,drop=F];
  x_curr_orig = x_all_orig[(n_hist + 1):(n_hist + n_curr),,drop=F];
  x_new_orig = x_all_orig[(n_hist + n_curr + 1):(n_hist + n_curr + n_new),,drop=F];
  #
  x_hist_aug = x_all_aug[1:n_hist,,drop=F];
  x_curr_aug = x_all_aug[(n_hist + 1):(n_hist + n_curr),,drop=F];
  x_new_aug = x_all_aug[(n_hist + n_curr + 1):(n_hist + n_curr + n_new),,drop=F];
  #
  y_hist = y_all[1:n_hist];
  y_curr = y_all[(n_hist+1):(n_hist + n_curr)];
  y_new = y_all[(n_hist + n_curr + 1):(n_hist + n_curr + n_new)];

  risk_new = risk_all[(n_hist + n_curr + 1):(n_hist + n_curr + n_new)];

  list(x_hist_orig = x_hist_orig,
       x_hist_aug = x_hist_aug,
       y_hist = y_hist,
       x_curr_orig = x_curr_orig,
       x_curr_aug = x_curr_aug,
       y_curr = y_curr,
       x_new_orig = x_new_orig,
       x_new_aug = x_new_aug,
       y_new = y_new,
       lin_pred_x_orig = lin_pred_x_orig,
       lin_pred_x_aug = lin_pred_x_aug,
       risk_new = risk_new)

}
