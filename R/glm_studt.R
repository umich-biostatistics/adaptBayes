#' Fit GLM with a regularized student-t prior regression coefficients
#'
#' Program for fitting a GLM equipped with a regularized student-t prior on the regression
#' coefficients, parametrized using the normal-inverse-gamma distribution. The 'regularization'
#' refers to the fact that the inverse-gamma scale is has a finite upper bound that
#' it smoothly approaches. This method was not used in the simulation study but was
#' used in the data analysis. Specifically, it corresponds to 'PedRESC2'.
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
#' @param beta_scale (pos. real) constants indicating the prior scale of the student-t prior.
#' @param dof (pos. integer) degrees of freedom for the student-t prior
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
#' @param return_as_stanfit (logical) should the function return the stanfit
#' object asis or should a summary of stanfit be returned as a regular list
#'
#' @return \code{list} object containing the draws and other information.
#'
#' @examples
#'
#'
#'     data(historical)
#'
#'     foo = glm_studt(y = historical$y_hist,
#'                     x_standardized = historical[,2:5],
#'                     family = "binomial",
#'                     beta_scale = 0.0231,
#'                     dof = 1,
#'                     slab_precision = 0.00444,
#'                     only_prior = 0,
#'                     mc_warmup = 1000,
#'                     mc_iter_after_warmup = 1000,
#'                     mc_chains = 2,
#'                     mc_thin = 1,
#'                     mc_stepsize = 0.1,
#'                     mc_adapt_delta = 0.99,
#'                     mc_max_treedepth = 15);
#'
#'     data(current)
#'
#'     foo = glm_studt(y = current$y_curr,
#'                     x_standardized = current[,2:11],
#'                     family = "binomial",
#'                     beta_scale = 0.0231,
#'                     dof = 1,
#'                     slab_precision = 0.00444,
#'                     only_prior = 0,
#'                     mc_warmup = 1000,
#'                     mc_iter_after_warmup = 1000,
#'                     mc_chains = 2,
#'                     mc_thin = 1,
#'                     mc_stepsize = 0.1,
#'                     mc_adapt_delta = 0.99,
#'                     mc_max_treedepth = 15);
#'
#' @import cmdstanr dplyr
#'
#' @export

glm_studt = function(y,
                     x_standardized,
                     family = "binomial",
                     beta_scale,
                     dof = 1,
                     slab_precision = (1/15)^2,
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

  # Now we do the sampling in Stan
  model_file <-
    system.file("stan",
                paste0("regstudt_", family, ".stan"),
                package = "adaptBayes",
                mustWork = TRUE)
  model <- cmdstanr::cmdstan_model(model_file)

  curr_fit <-
    tryCatch.W.E(
      model$sample(
        data = list(n_stan = length(y),
                    p_stan = ncol(x_standardized),
                    y_stan = y,
                    x_standardized_stan = x_standardized,
                    dof_stan = dof,
                    beta_scale_stan = beta_scale,
                    slab_precision_stan = slab_precision,
                    only_prior = as.integer(only_prior)),
        iter_warmup = mc_warmup,
        iter = mc_iter_after_warmup,
        chains = mc_chains,
        parallel_chains = min(mc_chains, getOption("mc.cores")),
        thin = mc_thin,
        step_size = mc_stepsize,
        adapt_delta = mc_adapt_delta,
        max_treedepth = mc_max_treedepth));

  if("simpleError"%in%class(curr_fit$value) || "error"%in%class(curr_fit$value)) {
    stop(curr_fit$value);
  }


  if(return_as_stanfit) {
    curr_fit$value;

  } else {

    model_diagnostics <- curr_fit$value$sampler_diagnostics()
    model_summary <- curr_fit$value$summary()

    list(num_divergences = sum(model_diagnostics[,,"divergent__"]),
         max_rhat = max(model_summary$rhat, na.rm=T),
         curr_beta0 = curr_fit$value$draws("mu", format="matrix")[, 1, drop = T],
         curr_beta = curr_fit$value$draws("beta", format="matrix"),
         theta = curr_fit$value$draws("theta", format="matrix"));
  }
}

