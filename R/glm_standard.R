#' GLM equipped with the 'standard' prior evaluated
#'
#' Program for fitting a GLM equipped with the 'standard' prior evaluated
#' in Boonstra and Barbaro, which is the regularized horseshoe.
#'
#' @param y (vector) outcomes corresponding to the type of glm desired. This should
#' match whatever datatype is expected by the stan program.
#' @param x_standardized (matrix) matrix of numeric values with number of rows equal
#' to the length of y and number of columns equal to p+q. It is assumed without
#' verification that each column is standardized to whatever scale the prior
#' expects - in Boonstra and Barbaro, all predictors are marginally generated to have
#' mean zero and unit variance, so no standardization is conducted. In practice,
#' all data should be standardized to have a common scale before model fitting.
#' If regression coefficients on the natural scale are desired, they be easily obtained
#' through unstandardizing.
#' @param family (character) Similar to argument in `glm` with the same name, but
#'  here this must be a character, and currently only 'binomial' (if y is binary) or
#' 'gaussian' (if y is continuous) are valid choices.
#' @param p see `q`
#' @param q (nonneg. integers) numbers, the sum of which add up to the number of columns
#' in x_standardized. For the standard prior, this distinction is only needed if a different
#' constant scale parameter (beta_orig_scale, beta_aug_scale), which is the constant 'c'
#' in the notation of Boonstra and Barbaro, is used.
#' @param beta_orig_scale see `beta_aug_scale`
#' @param beta_aug_scale (pos. real) constants indicating the prior scale of the
#' horseshoe. Both values correspond to 'c / sigma' in the notation of Boonstra and Barbaro,
#' because that paper never considers beta_orig_scale!=beta_aug_scale. Use
#' the function `solve_for_hiershrink_scale` to calculate this quantity. If 'y'
#' is binary, then sigma doesn't actually exist as a
#' parameter, and it will be set equal to 2 inside the function.
#' If 'y' is continuous, then sigma is equipped with its own weak
#' prior. In either case, it is not intended that the user scale by sigma "manually".
#' @param local_dof (pos. integer) number indicating the degrees of freedom for
#' lambda_j. Boonstra and Barbaro always used local_dof = 1. Choose a negative
#' value to tell the function that there are no local hyperparameters.
#' @param global_dof (pos. integer) number indicating the degrees of freedom for
#' tau. Boonstra and Barbaro always used global_dof = 1. Choose a negative
#' value to tell the function that there is no global hyperparameter.
#' @param slab_precision (pos. real) the slab-part of the regularized horseshoe,
#' this is equivalent to (1/d)^2 in the notation of Boonstra and Barbaro
#' @param intercept_offset (vector) vector of 0's and 1's equal having the same length as y.
#' Those observations with a value of 1 have an additional constant offset in their linear
#' predictor, effectively a different intercept. This is useful to jointly regress
#' two datasets in which it is believed that the regression coefficients are the same
#' but not the intercepts and could be useful (but was not used) in the simulation study
#' to compare to a benchmark, namely if both the historical and current datasets
#' were available but there is a desire to adjust for potentially different baseline prevalences.
#' @param only_prior (logical) should all data be ignored, sampling only from the prior?
#' @param mc_warmup number of MCMC warm-up iterations
#' @param mc_iter_after_warmup number of MCMC iterations after warm-up
#' @param mc_chains number of MCMC chains
#' @param mc_thin every nth draw to keep
#' @param mc_stepsize positive stepsize
#' @param mc_adapt_delta between 0 and 1
#' @param mc_max_treedepth max tree depth
#' @param return_as_stanfit (logical) should the function return the stanfit
#' object asis or should a summary of stanfit be returned as a regular list
#'
#' @import cmdstanr dplyr
#'
#' @return \code{list} object containing the draws and other information.
#'
#' @examples
#'
#' data(historical)
#'
#' foo = glm_standard(y = historical$y_hist,
#'                    x_standardized = historical[,2:5],
#'                    family = "binomial",
#'                    p = 4,
#'                    q = 0,
#'                    beta_orig_scale = 0.0231,
#'                    beta_aug_scale = 0.0231,
#'                    local_dof = 1,
#'                    global_dof = 1,
#'                    slab_precision = 0.00444,
#'                    intercept_offset = NULL,
#'                    only_prior = 0,
#'                    mc_warmup = 1000,
#'                    mc_iter_after_warmup = 1000,
#'                    mc_chains = 2,
#'                    mc_thin = 1,
#'                    mc_stepsize = 0.1,
#'                    mc_adapt_delta = 0.99,
#'                    mc_max_treedepth = 15);
#'
#'  data(current)
#'
#'  foo = glm_standard(y = current$y_curr,
#'                     x_standardized = current[,2:11],
#'                     family = "binomial",
#'                     p = 4,
#'                     q = 6,
#'                     beta_orig_scale = 0.0223,
#'                     beta_aug_scale = 0.0223,
#'                     local_dof = 1,
#'                     global_dof = 1,
#'                     slab_precision = 0.00444,
#'                     intercept_offset = NULL,
#'                     only_prior = 0,
#'                     mc_warmup = 1000,
#'                     mc_iter_after_warmup = 1000,
#'                     mc_chains = 2,
#'                     mc_thin = 1,
#'                     mc_stepsize = 0.1,
#'                     mc_adapt_delta = 0.99,
#'                     mc_max_treedepth = 15);
#'
#' @export

glm_standard = function(y,
                        x_standardized,
                        family = "binomial",
                        p,
                        q,
                        beta_orig_scale,
                        beta_aug_scale,
                        local_dof = 1,
                        global_dof = 1,
                        slab_precision = (1/15)^2,
                        intercept_offset = NULL,
                        only_prior = F,
                        mc_warmup = 50,
                        mc_iter_after_warmup = 50,
                        mc_chains = 1,
                        mc_thin = 1,
                        mc_stepsize = 0.1,
                        mc_adapt_delta = 0.9,
                        mc_max_treedepth = 15,
                        return_as_stanfit = FALSE) {

  if(family != "gaussian" && family != "binomial") {
    stop("'family' must equal 'gaussian' or 'binomial'")
  }

  stopifnot(ncol(x_standardized) == (p+q));
  if(is.null(intercept_offset)) {intercept_offset = numeric(length(y));}



  # Now we do the sampling in Stan
  model_file <-
    system.file("stan",
                paste0("reghs_", family, ".stan"),
                package = "adaptBayes",
                mustWork = TRUE)
  model <- cmdstanr::cmdstan_model(model_file)

  curr_fit <-
    tryCatch.W.E(
      model$sample(
        data = list(n_stan = length(y),
                    p_stan = p,
                    q_stan = q,
                    y_stan = y,
                    x_standardized_stan = x_standardized,
                    local_dof_stan = local_dof,
                    global_dof_stan = global_dof,
                    beta_orig_scale_stan = beta_orig_scale,
                    beta_aug_scale_stan = beta_aug_scale,
                    slab_precision_stan = slab_precision,
                    intercept_offset_stan = intercept_offset,
                    only_prior = as.integer(only_prior)),
        iter_warmup = mc_warmup,
        iter = mc_iter_after_warmup,
        chains = mc_chains,
        parallel_chains = min(mc_chains, getOption("mc.cores")),
        thin = mc_thin,
        step_size = mc_stepsize,
        adapt_delta = mc_adapt_delta,
        max_treedepth = mc_max_treedepth))

  if("simpleError"%in%class(curr_fit$value) || "error"%in%class(curr_fit$value)) {
    stop(curr_fit$value);
  }

  if(return_as_stanfit) {
    curr_fit$value;

  } else {
    model_diagnostics <- curr_fit$value$sampler_diagnostics()
    model_summary <- curr_fit$value$summary()

    if(q > 0) {

      list(num_divergences = sum(model_diagnostics[,,"divergent__"]),
           max_rhat = max(model_summary$rhat, na.rm=T),
           hist_beta0 = curr_fit$value$draws("mu", format="matrix")[, 1, drop = T],
           curr_beta0 = curr_fit$value$draws("mu", format="matrix")[, 1, drop = T] +
             curr_fit$value$draws("mu_offset", format="matrix")[, 1, drop = T],
           curr_beta = curr_fit$value$draws("beta", format="matrix"),
           theta_orig =  curr_fit$value$draws("theta_orig", format="matrix"),
           theta_aug = curr_fit$value$draws("theta_aug", format="matrix"));
    } else {
      list(num_divergences = sum(model_diagnostics[,,"divergent__"]),
           max_rhat = max(model_summary$rhat, na.rm=T),
           hist_beta0 = curr_fit$value$draws("mu", format="matrix")[, 1, drop = T],
           curr_beta0 = curr_fit$value$draws("mu", format="matrix")[, 1, drop = T] +
             curr_fit$value$draws("mu_offset", format="matrix")[, 1, drop = T],
           curr_beta = curr_fit$value$draws("beta", format="matrix"),
           theta_orig =  curr_fit$value$draws("theta_orig", format="matrix"));
    }
  }
}
