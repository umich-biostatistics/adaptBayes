#' Compute inverse of 'solve_for_hiershrink_scale'
#'
#' Instead of providing a desired effective number of parameters, the user provides
#' the scale value, which is c / sigma, in the notation of Boonstra and Barbaro, and the
#' function gives the implied prior number of effective parameters based upon this value.
#'
#' @param scale (pos. real) the scale parameter, defined as c / sigma.
#' @param npar (pos. integers): the number of covariates.
#' @param local_dof (pos. integer) number indicating the degrees of freedom for
#' lambda_j. Boonstra and Barbaro always used local_dof = 1. Choose a negative
#' value to tell the function that there are no local hyperparameters.
#' @param global_dof (pos. integer) number indicating the degrees of freedom for
#' tau. Boonstra and Barbaro always used global_dof = 1. Choose a negative
#' value to tell the function that there is no global hyperparameter.
#' @param slab_precision (pos. real) the slab-part of the regularized horseshoe,
#' this is equivalent to (1/d)^2 in the notation of Boonstra and Barbaro
#' @param n sample size
#' @param n_sim number of simulates
#'
#'
#' @return the implied number of effective parameters.
#'
#' @importFrom stats rt
#' @export

calculate_m_eff = function(scale,
                           npar,
                           local_dof = 1,
                           global_dof = 1,
                           slab_precision = (1/15)^2,
                           n,
                           n_sim = 2e5
) {
  do_local = (local_dof > 0);
  do_global = (global_dof > 0);
  stopifnot(isTRUE(all.equal(npar%%1,0)));#Ensure integers
  if(do_local) {
    lambda = matrix(rt(n_sim * npar,df = local_dof), nrow = n_sim);
  } else {
    lambda = matrix(1, nrow = n_sim, ncol = npar);
  }
  if(do_global) {
    tau = rt(n_sim,df=global_dof);
    lambda = lambda * tau;
    rm(tau);
  }
  random_scales = 1 / (slab_precision + 1/(scale^2 * lambda^2));
  kappa = 1/(1 + n * random_scales);
  prior_num = mean(rowSums(1 - kappa))

  return(prior_num);
}
