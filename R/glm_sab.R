#' Fit GLM with the 'sensible adaptive bayes' prior
#'
#' Program for fitting a GLM equipped with the 'sensible adaptive bayes' prior
#' evaluated in the manuscript.
#'
#' @param stan_fit an R object of class stanfit, which allows the function to run
#' without recompiling the stan code.
#' @param y (vector) outcomes corresponding to the type of glm desired. This should
#' match whatever datatype is expected by the stan program.
#' @param x_standardized (matrix) matrix of numeric values with number of rows equal
#' to the length of y and number of columns equal to p+q. It is assumed without
#' verification that each column is standardized to whatever scale the prior
#' expects - in Boonstra and Barbaro, all predictors are marginally generated to have
#' mean zero and unit variance, so no standardization is conducted. In practice,
#' all data should be standardized to have a common scale before model fitting.
#' If regression coefficients on the natural scale are desired, they can be easily obtained
#' through unstandardizing.
#' @param alpha_prior_mean (vector) p-length vector giving the mean of alpha from the
#' historical analysis, corresponds to m_alpha in Boonstra and Barbaro
#' @param alpha_prior_cov (matrix) pxp positive definite matrix giving the variance of
#' alpha from the historical analysis, corresponds to S_alpha in Boonstra and Barbaro
#' @param aug_projection (matrix) pxq matrix that approximately projects the regression coefficients of
#' the augmented predictors onto the space of the regression coefficients for the original
#' predictors.This is the matrix P in the notation of Boonstra and Barbaro.
#' It can be calculated using the function 'create_projection'
#' @param alpha_prior_cov (matrix) pxp positive definite matrix giving the variance of
#' alpha from the historical analysis, corresponds to S_alpha in Boonstra and Barbaro
#' @param phi_mean (real) mean of phi corresponding to a truncated normal distribution.
#' Since the support of the distribution is truncated to [0,1], it would make sense,
#' but is not required, that 'phi_mean' itself also be in [0,1]
#' @param beta_orig_scale,
#' @param beta_aug_scale (pos. real) constants indicating the prior scale of the
#' horseshoe. Both values correspond to 'c' in the notation of Boonstra and Barbaro,
#' because that paper never considers beta_orig_scale!=beta_aug_scale
#' @param local_dof (pos. integer) numbers indicating the degrees of freedom for
#' lambda_j and tau, respectively. Boonstra, et al. never considered local_dof != 1
#' or global_dof != 1.
#' @param global_dof (pos. integer) numbers indicating the degrees of freedom for
#' lambda_j and tau, respectively. Boonstra, et al. never considered local_dof != 1
#' or global_dof != 1.
#' @param slab_precision (pos. real) the slab-part of the regularized horseshoe,
#' this is equivalent to (1/d)^2 in the notation of Boonstra and Barbaro
#' @param only_prior (logical) should all data be ignored, sampling only from the prior?
#' @param mc_warmup number of MCMC warm-up iterations
#' @param mc_iter_after_warmup number of MCMC iterations after warm-up
#' @param mc_chains number of MCMC chains
#' @param mc_thin every nth draw to keep
#' @param mc_stepsize positive stepsize
#' @param mc_adapt_delta between 0 and 1
#' @param mc_max_treedepth max tree depth
#' @param ntries (pos. integer) the stan function will run up to this many times,
#' stopping either when the number of divergent transitions* is zero or when ntries
#' has been reached. The reported fit will be that with the fewest number of divergent iterations.
#' @param return_as_stanfit (logical) should the function return the stanfit
#' object asis or should a summary of stanfit be returned as a regular list
#' @param eigendecomp_hist_var: R object of class 'eigen' containing a pxp matrix
#' of eigenvectors in each row (equivalent to v_0 in Boonstra and Barbaro) and
#' a p-length vector of eigenvalues. This is by default equal to eigen(alpha_prior_cov)
#' @param scale_to_variance225: a vector assumed to be such that, when multiplied
#' by the diagonal elements of alpha_prior_cov, the result is a vector of
#' elements each equal to 225. This is explicitly calculated if it is not provided
#'
#' @return \code{list} object containing the draws and other information.
#'
#' @examples
#'
#' data(current)
#'
#' alpha_prior_cov = matrix(data = c(0.02936, -0.02078, 0.00216, -0.00637,
#'                                   -0.02078, 0.03192, -0.01029, 0.00500,
#'                                   0.00216, -0.01029, 0.01991, -0.00428,
#'                                   -0.00637, 0.00500, -0.00428, 0.01650),
#'                          byrow = FALSE, nrow = 4);
#'
#' scale_to_variance225 = diag(alpha_prior_cov) / 225;
#' eigendecomp_hist_var = eigen(alpha_prior_cov);
#' aug_projection1 = matrix(data = c(0.0608, -0.02628, -0.0488, 0.0484, 0.449, -0.0201,
#'                                   0.5695, -0.00855, 0.3877, 0.0729, 0.193, 0.4229,
#'                                   0.1816, 0.37240, 0.1107, 0.1081, -0.114, 0.3704,
#'                                   0.1209, 0.03683, -0.1517, 0.2178, 0.344, -0.1427),
#'                          byrow = TRUE, nrow = 4);
#'
#' foo = glm_sab(y = current$y_curr,
#'               x_standardized = current[,2:11],
#'               alpha_prior_mean = c(1.462, -1.660, 0.769, -0.756),
#'               alpha_prior_cov = alpha_prior_cov,
#'               aug_projection = aug_projection1,
#'               phi_mean = 1,
#'               phi_sd = 0.25,
#'               beta_orig_scale = 0.0223,
#'               beta_aug_scale = 0.0223,
#'               local_dof = 1,
#'               global_dof = 1,
#'               slab_precision = 0.00444,
#'               only_prior = 0,
#'               mc_warmup = 1000,
#'               mc_iter_after_warmup = 1000,
#'               mc_chains = 2,
#'               mc_thin = 1,
#'               mc_stepsize = 0.1,
#'               mc_adapt_delta = 0.999,
#'               mc_max_treedepth = 15,
#'               ntries = 2,
#'               eigendecomp_hist_var = eigendecomp_hist_var,
#'               scale_to_variance225 = scale_to_variance225);
#'
#' @import rstan
#' @export

glm_sab = function(stan_fit = stanmodels$SAB_Stable,
                   y,
                   x_standardized,
                   alpha_prior_mean,
                   alpha_prior_cov,
                   aug_projection,
                   phi_mean,
                   phi_sd,
                   beta_orig_scale,
                   beta_aug_scale,
                   local_dof = 1,
                   global_dof = 1,
                   slab_precision = (1/15)^2,
                   only_prior = F,
                   mc_warmup = 50,
                   mc_iter_after_warmup = 50,
                   mc_chains = 1,
                   mc_thin = 1,
                   mc_stepsize = 0.1,
                   mc_adapt_delta = 0.9,
                   mc_max_treedepth = 15,
                   ntries = 1,
                   return_as_stanfit = FALSE,
                   eigendecomp_hist_var = NULL,
                   scale_to_variance225 = NULL
) {

  if(is.null(eigendecomp_hist_var)) {
    eigendecomp_hist_var = eigen(alpha_prior_cov);
  }
  eigenvec_hist_var = t(eigendecomp_hist_var$vectors);
  sqrt_eigenval_hist_var = sqrt(eigendecomp_hist_var$values);

  if(is.null(scale_to_variance225)) {
    scale_to_variance225 = diag(alpha_prior_cov) / 225;
  }

  p = length(alpha_prior_mean);
  q = ncol(x_standardized) - p;
  if(p == 1) {
    alpha_prior_mean = array(alpha_prior_mean,dim=1);
    sqrt_eigenval_hist_var = array(sqrt_eigenval_hist_var,dim=1);
    scale_to_variance225 = array(scale_to_variance225,dim=1);
  }

  max_divergences = -Inf;
  accepted_divergences = Inf;
  curr_try = 1;

  while(curr_try <= ntries) {
    assign("curr_fit",tryCatch.W.E(sampling(object = stan_fit,
                                            data = list(n_stan = length(y),
                                                        p_stan = p,
                                                        q_stan = q,
                                                        y_stan = y,
                                                        x_standardized_stan = x_standardized,
                                                        aug_projection_stan = aug_projection,
                                                        alpha_prior_mean_stan = alpha_prior_mean,
                                                        alpha_prior_cov_stan = alpha_prior_cov,
                                                        sqrt_eigenval_hist_var_stan = sqrt_eigenval_hist_var,
                                                        eigenvec_hist_var_stan = eigenvec_hist_var,
                                                        local_dof_stan = local_dof,
                                                        global_dof_stan = global_dof,
                                                        beta_orig_scale_stan = beta_orig_scale,
                                                        beta_aug_scale_stan = beta_aug_scale,
                                                        slab_precision_stan = slab_precision,
                                                        scale_to_variance225 = scale_to_variance225,
                                                        phi_mean_stan = phi_mean,
                                                        phi_sd_stan = phi_sd,
                                                        only_prior = as.integer(only_prior)),
                                           warmup = mc_warmup,
                                           iter = mc_iter_after_warmup + mc_warmup,
                                           chains = mc_chains,
                                           thin = mc_thin,
                                           control = list(stepsize = mc_stepsize,
                                                          adapt_delta = mc_adapt_delta,
                                                          max_treedepth = mc_max_treedepth))));

    if("simpleError"%in%class(curr_fit$value) || "error"%in%class(curr_fit$value)) {
      stop(curr_fit$value);
    }
    if(return_as_stanfit) {
      break;
    }
    curr_divergences = count_stan_divergences(curr_fit$value);
    rhat_check = max(summary(curr_fit$value)$summary[,"Rhat"],na.rm=T);
    # Originally, the break conditions were baesd upon having both no divergent
    # transitions as well as a max Rhat (i.e. gelman-rubin diagnostic) sufficiently
    # close to 1. I subsequently changed the conditions to be based only upon the
    # first, which is reflected by setting rhat = T immediately below.
    break_conditions = c(divergence = F, rhat = T);
    if(curr_divergences == 0) {
      max_divergences = 0;
      break_conditions["divergence"] = T;
    } else {
      max_divergences = max(max_divergences, curr_divergences, na.rm = T);
      curr_try = curr_try + 1;
    }
    #update if fewer divergent transitions were found
    if(curr_divergences < accepted_divergences) {
      accepted_divergences = curr_divergences;
      max_rhat = rhat_check;
      foo = rstan::extract(curr_fit$value);
      curr_beta0 = as.numeric(foo$mu);
      curr_beta = foo$beta;
      theta_orig = foo$theta_orig;
      theta_aug = foo$theta_aug;
      phi = foo$phi_copy;
      eta = foo$eta;
    }
    if(all(break_conditions)) {
      break;
    }
  }
  if(return_as_stanfit) {
    curr_fit$value;
  } else {
    list(accepted_divergences = accepted_divergences,
         max_divergences = max_divergences,
         max_rhat = max_rhat,
         curr_beta0 = curr_beta0,
         curr_beta = curr_beta,
         theta_orig = theta_orig,
         theta_aug = theta_aug,
         phi = phi,
         eta = eta);
  }
}
