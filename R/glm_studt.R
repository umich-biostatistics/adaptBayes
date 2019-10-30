#' Fit GLM with a regularized student-t prior regression coefficients
#'
#' Program for fitting a GLM equipped with a regularized student-t prior on the regression
#' coefficients, parametrized using the normal-inverse-gamma distribution. The 'regularization'
#' refers to the fact that the inverse-gamma scale is has a finite upper bound that
#' it smoothly approaches. This method was not used in the simulation study but was
#' used in the data analysis. Specifically, it corresponds to 'PedRESC2'.
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
#' If regression coefficients on the natural scale are desired, they be easily obtained
#' through unstandardizing.
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
#' @param ntries (pos. integer) the stan function will run up to this many times,
#' stopping either when the number of divergent transitions* is zero or when ntries
#' has been reached. The reported fit will be that with the fewest number of divergent iterations.
#' @param return_as_stanfit (logical) should the function return the stanfit
#' object asis or should a summary of stanfit be returned as a regular list
#'
#' @return \code{list} object containing the draws and other information.
#'
#' @import rstan
#'
#' @export

glm_studt = function(stan_fit = stanmodels$RegStudT,
                     y,
                     x_standardized,
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
                     ntries = 1,
                     return_as_stanfit = FALSE) {

  max_divergences = -Inf;
  accepted_divergences = Inf;
  curr_try = 1;

  while(curr_try <= ntries) {
    assign("curr_fit",tryCatch.W.E(sampling(object = stan_fit,
                                           data = list(n_stan = length(y),
                                                       p_stan = ncol(x_standardized),
                                                       y_stan = y,
                                                       x_standardized_stan = x_standardized,
                                                       dof_stan = dof,
                                                       beta_scale_stan = beta_scale,
                                                       slab_precision_stan = slab_precision,
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
      theta = foo$theta;
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
         theta = theta);
  }
}
