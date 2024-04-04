#' Fit GLM with a regularized student-t prior regression coefficients
#'
#' Program for fitting a GLM equipped with a regularized student-t prior on the
#' regression coefficients, parametrized using the normal-inverse-gamma
#' distribution. The 'regularization' refers to the fact that the inverse-gamma
#' scale is has a finite upper bound that it smoothly approaches. This method
#' was not used in the simulation study but was used in the data analysis.
#' Specifically, it corresponds to 'PedRESC2'.
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
#'   are desired, they be easily obtained through unstandardizing.
#' @param family (character) Similar to argument in `glm` with the same name,
#'   but here this must be a character, and currently only 'binomial' (if y is
#'   binary) or 'gaussian' (if y is continuous) are valid choices.
#' @param beta_scale (pos. real) constants indicating the prior scale of the
#'   student-t prior.
#' @param dof (pos. integer) degrees of freedom for the student-t prior
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
#'                     mu_sd = 5,
#'                     only_prior = 0,
#'                     mc_warmup = 200,
#'                     mc_iter_after_warmup = 200,
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
#'                     mu_sd = 5,
#'                     only_prior = 0,
#'                     mc_warmup = 200,
#'                     mc_iter_after_warmup = 200,
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
                     seed = sample.int(.Machine$integer.max, 1),
                     slab_precision = NULL
) {


  if(family != "gaussian" && family != "binomial") {
    stop("'family' must equal 'gaussian' or 'binomial'")
  }

  if(!is.null(slab_precision)) {
    slab_dof = Inf;
    slab_scale = 1 / sqrt(slab_precision);
    message(paste0("'slab_precision' will be going away in a future release; use 'slab_dof = Inf' and 'slab_scale = 1/sqrt(",slab_precision,")'"))
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
                    slab_dof_stan = slab_dof,
                    slab_scale_stan = slab_scale,
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
        refresh = 0));

  if("simpleError"%in%class(curr_fit$value) || "error"%in%class(curr_fit$value)) {
    stop(curr_fit$value);
  }


  if(return_as_CmdStanMCMC) {
    curr_fit$value;

  } else {

    model_diagnostics <- curr_fit$value$diagnostic_summary()
    model_summary <- curr_fit$value$summary()

    list(num_divergences = sum(model_diagnostics$num_divergent),
         num_max_treedepth = sum(model_diagnostics$num_max_treedepth),
         min_ebfmi = min(model_diagnostics$ebfmi),
         max_rhat = max(model_summary$rhat, na.rm=T),
         mu = curr_fit$value$draws("mu", format="matrix")[, 1, drop = T],
         beta = curr_fit$value$draws("beta", format="matrix"),
         theta = curr_fit$value$draws("theta", format="matrix"),
         slab = curr_fit$value$draws("slab_copy", format="matrix"));
  }
}

