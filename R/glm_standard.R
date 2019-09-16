
#DESCRIPTION: Program for fitting a GLM equipped with the 'standard' prior evaluated
#in Boonstra and Barbaro, which is the regularized horseshoe. It has two intended uses:
#compile stan scripts or analyze data. First, if the user provides nothing but a valid
#stan_path, then the stan script is compiled. Second, the user provides both a compiled
#stanfit object as well asvalues for y, x_standardized, #, q, and any other desired
#arguments to actually fit a regression.
#
#
#ARGUMENTS:
#
#stan_fit: an R object of class stanfit, which allows the function to run without recompiling the stan code.
#
#stan_path: (character) a path pointing to a .stan file, which indicates the stan code to compile and run. If
#both stan_fit and stan_path are provided, stan_fit takes precedence.
#
#y (vector) outcomes corresponding to the type of glm desired. This should match whatever datatype is expected
#by the stan program.
#
#x_standardized (matrix) matrix of numeric values with number of rows equal to the length of y and number of columns
#equal to p+q. It is assumed without verification that each column is standardized to whatever scale the prior
#expects - in Boonstra and Barbaro, all predictors are marginally generated to have mean zero and unit variance, so no
#standardization is conducted. In practice, all data should be standardized to have a common scale before model fitting.
#If regression coefficients on the natural scale are desired, they be easily obtained through unstandardizing.
#
#p, q (nonneg. integers) numbers, the sum of which add up to the number of columns in x_standardized. For the standard
#prior, this distinction is only needed if a different constant scale parameter (beta_orig_scale, beta_aug_scale), which is
#the constant 'c' in the notation of Boonstra and Barbaro, is used.
#
#beta_orig_scale, beta_aug_scale (pos. real) constants indicating the prior scale of the horseshoe. Both values correspond
#to 'c' in the notation of Boonstra and Barbaro, because that paper never considers beta_orig_scale!=beta_aug_scale
#
#local_dof, global_dof (pos. integer) numbers indicating the degrees of freedom for lambda_j and tau, respectively. Boonstra,
#et al. never considered local_dof != 1 or global_dof != 1.
#
#slab_precision (pos. real) the slab-part of the regularized horseshoe, this is equivalent to (1/d)^2 in the notation of
#Boonstra and Barbaro
#
#intercept_offset (vector) vector of 0's and 1's equal having the same length as y. Those observations with a value of 1
#have an additional constant offset in their linear predictor, effectively a different intercept. This is useful to jointly
#regress two datasets in which it is believed that the regression coefficients are the same but not the intercepts and could be
#useful (but was not used) in the simulation study to compare to a benchmark, namely if both the historical and current datasets
#were available but there is a desire to adjust for potentially different baseline prevalences.
#
#only_prior (logical) should all data be ignored, sampling only from the prior?
#
#ntries (pos. integer) the stan function will run up to this many times, stopping either when the number of
#*divergent transitions* is zero or when ntries has been reached. The reported fit will be that with the fewest number of
#divergent iterations.

#' GLM equipped with the 'standard' prior evaluated
#'
#' Program for fitting a GLM equipped with the 'standard' prior evaluated
#' in Boonstra and Barbaro, which is the regularized horseshoe.
#' It has two intended uses: compile stan scripts or analyze data. First, if the
#' user provides nothing but a valid 'stan_path', then the stan script is compiled.
#' Second, the user provides both a compiled stanfit object as well asvalues for y,
#' x_standardized, #, q, and any other desired arguments to actually fit a regression.
#'
#' @param stan_fit an R object of class stanfit, which allows the function to run
#' without recompiling the stan code.
#' @param stan_path (character) a path pointing to a .stan file, which indicates
#' the stan code to compile and run. If both stan_fit and stan_path are provided,
#' stan_fit takes precedence.
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
#' @param p (nonneg. integers) numbers, the sum of which add up to the number of columns
#' in x_standardized. For the standard prior, this distinction is only needed if a different
#' constant scale parameter (beta_orig_scale, beta_aug_scale), which is the constant 'c'
#' in the notation of Boonstra and Barbaro, is used.
#' @param q (nonneg. integers) numbers, the sum of which add up to the number of columns
#' in x_standardized. For the standard prior, this distinction is only needed if a different
#' constant scale parameter (beta_orig_scale, beta_aug_scale), which is the constant 'c'
#' in the notation of Boonstra and Barbaro, is used.
#' @param beta_orig_scale (pos. real) constants indicating the prior scale of the
#' horseshoe. Both values correspond to 'c' in the notation of Boonstra and Barbaro,
#' because that paper never considers beta_orig_scale!=beta_aug_scale
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
#' @param ntries (pos. integer) the stan function will run up to this many times,
#' stopping either when the number of divergent transitions* is zero or when ntries
#' has been reached. The reported fit will be that with the fewest number of divergent iterations.
#'
#' @import rstan
#'
#' @return \code{list} object containing the draws and other information.
#'
#' @export

glm_standard = function(stan_fit = stanmodels$RegHS_stable,
                        stan_path,
                        y = c(0,1),
                        x_standardized = matrix(0,length(y),6),
                        p = 3,
                        q = 3,
                        beta_orig_scale = 1,
                        beta_aug_scale = 1,
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
                        ntries = 1) {

  stopifnot(ncol(x_standardized) == (p+q));

  max_divergences = -Inf;
  accepted_divergences = Inf;
  curr_try = 1;
  if(is.null(intercept_offset)) {intercept_offset = numeric(length(y));}

  while(curr_try <= ntries) {
    assign("curr_fit",tryCatch.W.E(stan(file = stan_path,
                                        fit = stan_fit,# ,
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
    if(!"stanfit"%in%class(stan_fit)) {
      break;
    }
    divergent_check = unlist(lapply(curr_fit$warning,grep,pattern="divergent transitions",value=T));
    rhat_check = max(summary(curr_fit$value)$summary[,"Rhat"],na.rm=T);
    #Originally, the break conditions were baesd upon having both no divergent transitions as well as a max Rhat (i.e. gelman-rubin
    #diagnostic) sufficiently close to 1. I subsequently changed the conditions to be based only upon the first, which is reflected
    #by setting rhat = T immediately below.
    break_conditions = c(divergence = F, rhat = T);
    if(length(divergent_check) == 0) {#corresponds to zero divergent transitions
      curr_divergences = 0;
      max_divergences = max(max_divergences,curr_divergences,na.rm=T);
      break_conditions["divergence"] = T;
    } else {#corresponds to > zero divergent transitions
      curr_divergences <- max(as.numeric(strsplit(divergent_check," ")$message),na.rm=T);
      max_divergences = max(max_divergences,curr_divergences,na.rm=T);
      curr_try = curr_try + 1;
    }
    #update if fewer divergent transitions were found
    if(curr_divergences < accepted_divergences) {
      accepted_divergences = curr_divergences;
      max_rhat = rhat_check;
      foo = rstan::extract(curr_fit$value);
      hist_beta0 = as.numeric(foo$mu);
      curr_beta0 = as.numeric(foo$mu) + as.numeric(foo$mu_offset);
      curr_beta = foo$beta;
      theta_orig = foo$theta_orig;
      theta_aug = foo$theta_aug;
    }
    if(all(break_conditions)) {
      break;
    }
  }
  if(!"stanfit"%in%class(stan_fit)) {
    curr_fit$value;
  } else {
    list(accepted_divergences = accepted_divergences,
         max_divergences = max_divergences,
         max_rhat = max_rhat,
         hist_beta0 = hist_beta0,
         curr_beta0 = curr_beta0,
         curr_beta = curr_beta,
         theta_orig = theta_orig,
         theta_aug = theta_aug);
  }
}
