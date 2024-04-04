#' Fit GLM with the 'naive adaptive bayes' prior
#'
#' Program for fitting a GLM equipped with the 'naive adaptive bayes' prior
#' evaluated in the manuscript.
#'
#' @param y (vector) outcomes corresponding to the type of glm desired. This
#'   should match whatever datatype is expected by the stan program.
#' @param x_standardized (matrix) matrix of numeric values with number of rows
#'   equal to the length of y and number of columns equal to p+q. It is assumed
#'   without verification that each column is standardized to whatever scale the
#'   prior expects - in Boonstra and Barbaro, all predictors are marginally
#'   generated to have mean zero and unit variance, so no standardization is
#'   conducted. In practice, all data should be standardized to have a common
#'   scale before model fitting. If regression coefficients on the natural scale
#'   are desired, they can be easily obtained through unstandardizing.
#' @param family (character) Similar to argument in `glm` with the same name,
#'   but here this must be a character, and currently only 'binomial' (if y is
#'   binary) or 'gaussian' (if y is continuous) are valid choices.
#' @param alpha_prior_mean (vector) p-length vector giving the mean of alpha
#'   from the historical analysis, corresponds to m_alpha in Boonstra and
#'   Barbaro
#' @param alpha_prior_cov (matrix) pxp positive definite matrix giving the
#'   variance of alpha from the historical analysis, corresponds to S_alpha in
#'   Boonstra and Barbaro
#' @param phi_dist (character) the name of the distribution to use as a prior on
#'   phi. This must be either 'trunc_norm' or 'beta'.
#' @param phi_mean see `phi_sd`
#' @param phi_sd (real) prior mean and standard deviation of phi. At a minimum,
#'   phi_mean must be between 0 and 1 (inclusive) and phi_sd must be
#'   non-negative (you *can* choose phi_sd = 0, meaning that phi is identically
#'   equal to phi_mean). If 'phi_dist' is 'trunc_norm', then 'phi_mean' and
#'   'phi_sd' are interpreted as the parameters of the *untruncated* normal
#'   distribution and so are not actually the parameters of the resulting
#'   distribution after truncating phi to the 0,1 interval. If 'phi_dist' is
#'   'beta', then 'phi_mean' and 'phi_sd' are interpreted as the literal mean
#'   and standard deviation, from which the shape parameters are calculated.
#'   When 'phi_dist' is 'beta', not all choices of 'phi_mean' and 'phi_sd' are
#'   valid, e.g. the standard deviation of the beta distribution must be no
#'   greater than sqrt(phi_mean * (1 - phi_mean)). Also, the beta distribution
#'   is difficult to sample from if one or both of the shape parameters is much
#'   less than 1. An error will be thrown if an invalid parameterization is
#'   provided, and a warning will be thrown if a parameterization is provided
#'   that is likely to result in a "challenging" prior.
#' @param eta_param (real) prior hyperparmeter for eta, which scales the
#'   `alpha_prior_cov` in the adaptive prior contribution and is apriori
#'   distributed as an inverse-gamma random variable. Specifically, `eta_param`
#'   is a common value for the shape and rate of the inverse-gamma, meaning that
#'   larger values cause the prior distribution of eta to concentrate around
#'   one. You may choose `eta_param = Inf` to make eta identically equal to 1
#' @param beta_orig_scale see `beta_aug_scale`
#' @param beta_aug_scale (pos. real) constants indicating the prior scale of the
#'   horseshoe. Both values correspond to 'c / sigma' in the notation of
#'   Boonstra and Barbaro, because that paper never considers
#'   beta_orig_scale!=beta_aug_scale. Use the function
#'   `solve_for_hiershrink_scale` to calculate this quantity. If 'y' is binary,
#'   then sigma doesn't actually exist as a parameter, and it will be set equal
#'   to 2 inside the function. If 'y' is continuous, then sigma is equipped with
#'   its own weak prior. In either case, it is not intended that the user scale
#'   by sigma "manually".
#' @param beta_aug_scale_tilde (pos. real) constant indicating the prior scale
#'   of the horseshoe for the augmented covariates when phi = 1, i.e. when the
#'   historical analysis is fully used. This corresponds to tilde_c in Boonstra
#'   and Barbaro
#' @param local_dof (pos. integer) number indicating the degrees of freedom for
#'   lambda_j. Boonstra and Barbaro always used local_dof = 1. Choose a negative
#'   value to tell the function that there are no local hyperparameters.
#' @param global_dof (pos. integer) number indicating the degrees of freedom for
#'   tau. Boonstra and Barbaro always used global_dof = 1. Choose a negative
#'   value to tell the function that there is no global hyperparameter.
#' @param slab_dof see `slab_scale`
#' @param slab_scale (pos. real) these control the slab-part of the regularized
#'   horseshoe. Specifically, in the notation of Boonstra and Barbaro,
#'   d^2~InverseGamma(`slab_dof`/2, `slab_scale`^2*`slab_dof`/2). In Boonstra and
#'   Barbaro, d was fixed at 15, and you can achieve this by leaving these at
#'   their default values of `slab_dof` = Inf and `slab_scale` = 15.
#' @param mu_sd (pos. real) the prior standard deviation for the intercept
#'   parameter mu
#' @param only_prior (logical) should all data be ignored, sampling only from
#'   the prior?
#' @param mc_warmup number of MCMC warm-up iterations
#' @param mc_iter_after_warmup number of MCMC iterations after warm-up
#' @param mc_chains number of MCMC chains
#' @param mc_thin every nth draw to keep
#' @param mc_stepsize positive stepsize
#' @param mc_adapt_delta between 0 and 1
#' @param mc_max_treedepth max tree depth
#' @param return_as_CmdStanMCMC (logical) should the function return the CmdStanMCMC
#'   object asis or should a summary of CmdStanMCMC be returned as a regular list
#' @param eigendecomp_hist_var R object of class 'eigen' containing a pxp matrix
#'   of eigenvectors in each row (equivalent to v_0 in Boonstra and Barbaro) and
#'   a p-length vector of eigenvalues. This is by default equal to
#'   eigen(alpha_prior_cov)
#' @param scale_to_variance225 a vector assumed to be such that, when multiplied
#'   by the diagonal elements of alpha_prior_cov, the result is a vector of
#'   elements each equal to 225. This is explicitly calculated if it is not
#'   provided
#' @param seed seed for the underlying STAN model to allow for reproducibility
#' @param slab_precision (pos. real) the slab-part of the regularized horseshoe,
#'   this is equivalent to (1/d)^2 in the notation of Boonstra and Barbaro. If
#'   specified, it is assumed that you want a fixed slab component and will take
#'   precedence over any provided values of `slab_dof` and `slab_scale`;
#'   `slab_precision` is provided for backwards compatibility but will be going
#'   away in a future release, and the proper way to specify a fixed slab
#'   component with with precision 1/d^2 for some number d is through `slab_dof
#'   = Inf` and `slab_scale = d`.
#'
#' @return `list` object containing the draws and other information.
#'
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
#' scale_to_variance225 = diag(alpha_prior_cov) / 225;
#' eigendecomp_hist_var = eigen(alpha_prior_cov);
#'
#' foo = glm_nab(y = current$y_curr,
#'               x_standardized = current[,2:11],
#'               family = "binomial",
#'               alpha_prior_mean = c(1.462, -1.660, 0.769, -0.756),
#'               alpha_prior_cov = alpha_prior_cov,
#'               phi_dist = "trunc_norm",
#'               phi_mean = 1,
#'               phi_sd = 0.25,
#'               eta_param = 2.5,
#'               beta_orig_scale = 0.0223,
#'               beta_aug_scale = 0.0223,
#'               beta_aug_scale_tilde = 0.05,
#'               local_dof = 1,
#'               global_dof = 1,
#'               mu_sd = 5,
#'               only_prior = 0,
#'               mc_warmup = 200,
#'               mc_iter_after_warmup = 200,
#'               mc_chains = 2,
#'               mc_thin = 1,
#'               mc_stepsize = 0.1,
#'               mc_adapt_delta = 0.999,
#'               mc_max_treedepth = 15,
#'               eigendecomp_hist_var = eigendecomp_hist_var,
#'               scale_to_variance225 = scale_to_variance225);
#'
#' @import cmdstanr dplyr
#'
#' @export

glm_nab = function(y,
                   x_standardized,
                   family = "binomial",
                   alpha_prior_mean,
                   alpha_prior_cov,
                   phi_dist = "trunc_norm",
                   phi_mean = 1,
                   phi_sd = 0.25,
                   eta_param = 2.5,
                   beta_orig_scale,
                   beta_aug_scale,
                   beta_aug_scale_tilde,
                   local_dof = 1,
                   global_dof = 1,
                   slab_dof = Inf,
                   slab_scale = 15,
                   mu_sd = 5,
                   only_prior = F,
                   mc_warmup = 1e3,
                   mc_iter_after_warmup = 1e3,
                   mc_chains = 1,
                   mc_thin = 1,
                   mc_stepsize = 0.1,
                   mc_adapt_delta = 0.9,
                   mc_max_treedepth = 15,
                   return_as_CmdStanMCMC = FALSE,
                   eigendecomp_hist_var = NULL,
                   scale_to_variance225 = NULL,
                   seed = sample.int(.Machine$integer.max, 1),
                   slab_precision = NULL
) {

  if(family != "gaussian" && family != "binomial") {
    stop("'family' must equal 'gaussian' or 'binomial'")
  }

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

  if(phi_mean < 0 || phi_mean > 1) {stop("'phi_mean' should be between 0 and 1")}
  if(phi_sd < 0) {stop("'phi_sd' must be non-negative")}
  if(eta_param < 0) {stop("'eta_param' must be non-negative")}
  if(mu_sd < 0) {stop("'mu_sd' must be non-negative")}

  if(!is.null(slab_precision)) {
    slab_dof = Inf;
    slab_scale = 1 / sqrt(slab_precision);
    message(paste0("'slab_precision' will be going away in a future release; use 'slab_dof = Inf' and 'slab_scale = 1/sqrt(",slab_precision,")'"))
  }

  if(phi_dist == "beta") {

    sd_mean_fraction =
      case_when(
        # If phi_sd = 0, then phi equals phi_mean always (even for phi_mean = 1)
        phi_sd == 0 ~ 0,
        TRUE ~ phi_sd / sqrt(phi_mean * (1 - phi_mean)))
    if(sd_mean_fraction >= 1) {
      stop("'phi_sd' must be less than 'sqrt(phi_mean*(1-phi_mean))' to yield a valid beta distribution.")
    } else if (sd_mean_fraction >= 0.85) {
      warning("'phi_sd' may be too large relative to 'sqrt(phi_mean*(1-phi_mean))' to provide stable inference. Consider decreasing 'phi_sd' or moving 'phi_mean' closer to 0.5")
    }

  } else if(phi_dist != "trunc_norm") {
    stop("'phi_dist' must equal 'trunc_norm' or 'beta'")
  }

  # Now we do the sampling in Stan
  if(phi_mean == 1 && phi_sd == 0 && is.infinite(eta_param)) {
    model_file <-
      system.file("stan",
                  paste0("nab_simple_", family, ".stan"),
                  package = "adaptBayes",
                  mustWork = TRUE)
  } else {
    model_file <-
      system.file("stan",
                  paste0("nab_", family, ".stan"),
                  package = "adaptBayes",
                  mustWork = TRUE)
  }

  model <- cmdstanr::cmdstan_model(model_file)


  curr_fit <-
    tryCatch.W.E(
      model$sample(
        data = list(n_stan = length(y),
                    p_stan = p,
                    q_stan = q,
                    y_stan = y,
                    x_standardized_stan = x_standardized,
                    alpha_prior_mean_stan = alpha_prior_mean,
                    alpha_prior_cov_stan = alpha_prior_cov,
                    sqrt_eigenval_hist_var_stan = sqrt_eigenval_hist_var,
                    eigenvec_hist_var_stan = eigenvec_hist_var,
                    local_dof_stan = local_dof,
                    global_dof_stan = global_dof,
                    beta_orig_scale_stan = beta_orig_scale,
                    beta_aug_scale_stan = beta_aug_scale,
                    beta_aug_scale_tilde_stan = beta_aug_scale_tilde,
                    slab_dof_stan = slab_dof,
                    slab_scale_stan = slab_scale,
                    scale_to_variance225 = scale_to_variance225,
                    phi_prior_type = ifelse(phi_dist == "trunc_norm", 1L, 0L),
                    phi_mean_stan = phi_mean,
                    phi_sd_stan = phi_sd,
                    eta_param_stan = eta_param,
                    mu_sd_stan = mu_sd,
                    only_prior = as.integer(only_prior)),
        iter_warmup = mc_warmup,
        iter = mc_iter_after_warmup,
        chains = mc_chains,
        parallel_chains = min(mc_chains, getOption("mc.cores")),
        thin = mc_thin,
        step_size = mc_stepsize,
        adapt_delta = mc_adapt_delta,
        max_treedepth = mc_max_treedepth,
        seed = seed,
        refresh = 0))

  if("simpleError"%in%class(curr_fit$value) || "error"%in%class(curr_fit$value)) {
    stop(curr_fit$value);
  }

  if(return_as_CmdStanMCMC) {
    curr_fit$value;

  } else {

    model_diagnostics <- curr_fit$value$diagnostic_summary()
    model_summary <- curr_fit$value$summary()

    if(phi_mean == 1 && phi_sd == 0 && is.infinite(eta_param)) {
      phi = rep(1, mc_iter_after_warmup * mc_chains);
      eta = rep(1, mc_iter_after_warmup * mc_chains);
    } else {
      phi = curr_fit$value$draws("phi_copy", format="matrix")[, 1, drop = T];
      eta = curr_fit$value$draws("eta_copy", format="matrix")[, 1, drop = T];
    }

    list(num_divergences = sum(model_diagnostics$num_divergent),
         num_max_treedepth = sum(model_diagnostics$num_max_treedepth),
         min_ebfmi = min(model_diagnostics$ebfmi),
         max_rhat = max(model_summary$rhat, na.rm=T),
         mu = curr_fit$value$draws("mu", format="matrix")[, 1, drop = T],
         beta = curr_fit$value$draws("beta", format="matrix"),
         theta_orig =  curr_fit$value$draws("theta_orig", format="matrix"),
         theta_aug = curr_fit$value$draws("theta_aug", format="matrix"),
         slab = curr_fit$value$draws("slab_copy", format="matrix"),
         phi = phi,
         eta = eta);
  }
}
