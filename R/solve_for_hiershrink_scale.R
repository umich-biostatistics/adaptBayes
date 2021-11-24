#' Numerical-based solution to the scale parameter
#'
#' This function calculates a numerical-based solution to the quantity c / sigma in
#' the equation at the end of first paragraph on page e50. It is intended to
#' be provided as the value for beta_orig_scale and beta_aug_scale in the functions
#' `glm_sab`, `glm_nab`, and `glm_standard`
#'
#' If the outcome of interest is binary, then sigma doesn't actually exist as a real
#' parameter, and it will be set equal to 2 inside `glm_sab`,  `glm_nab`, or `glm_standard`.
#' If the outcome of interest is continuous, then sigma is equipped with its own weak
#' prior. In either case, it is not intended that the user scale by sigma "manually".
#'
#'
#' @param target_mean (pos. reals): the desired prior number of effective parameters
#' (tilde xi_eff in Boonstra and Barbaro). An error will be thrown if target_mean > npar
#' @param npar (pos. integers): the number of covariates.
#' @param local_dof (pos. integer) number indicating the degrees of freedom for
#' lambda_j. Boonstra and Barbaro always used local_dof = 1. Choose a negative
#' value to tell the function that there are no local hyperparameters.
#' @param global_dof (pos. integer) number indicating the degrees of freedom for
#' tau. Boonstra and Barbaro always used global_dof = 1. Choose a negative
#' value to tell the function that there is no global hyperparameter.
#' @param slab_precision (pos. real) the slab-part of the regularized horseshoe,
#' this is equivalent to (1/d)^2 in the notation of Boonstra and Barbaro
#' @param n (pos. integer) sample size of the study
#' @param tol (pos. real) numerical tolerance for convergence of solution
#' @param max_iter (pos. integer) maximum number of iterations to run without
#' convergence before giving up
#' @param n_sim (pos. integer) number of simulated draws from the underlying student-t
#' hyperpriors to calculate the Monte Carlo-based approximation of the expectation.
#'
#'
#' @return A \code{list} containing the following named elements:
#' \itemize{
#'   \item{scale1}{the value of }
#'   \item{diff_from_target1}{}
#'   \item{iter1}{}
#'   \item{prior_num1}{}
#' }
#'
#'
#' @export

solve_for_hiershrink_scale = function(target_mean,
                                      npar,
                                      local_dof = 1,
                                      global_dof = 1,
                                      slab_precision = (1/15)^2,
                                      n,
                                      tol = .Machine$double.eps^0.5,
                                      max_iter = 100,
                                      n_sim = 2e5
) {

  stopifnot(isTRUE(all.equal(npar%%1,0)));#Ensure integers
  do_local = (local_dof > 0);
  do_global = (global_dof > 0);
  if(do_local) {
    lambda = matrix(rt(n_sim * npar,df = local_dof), nrow = n_sim);
  } else {
    lambda = matrix(1, nrow = n_sim, ncol = npar);
  }
  if(do_global) {
    tau = rt(n_sim,df=global_dof);
    lambda = lambda* tau;
    rm(tau);
  }
  lambda_sq = lambda^2

  stopifnot(target_mean > 0 && target_mean < npar);#Ensure proper bounds
  log_scale =
    diff_target =
    numeric(max_iter);
  log_scale[1] = log(target_mean/(npar - target_mean) / sqrt(n));
  random_scales = 1 / (slab_precision + 1/(exp(2*log_scale[1]) * lambda_sq));
  kappa = 1 / (1 + n * random_scales);
  diff_target[1] = mean(rowSums(1 - kappa)) - target_mean;
  log_scale[2] = 0.02 + log_scale[1];
  random_scales = 1 / (slab_precision + 1/(exp(2*log_scale[2]) * lambda_sq));
  kappa = 1 / (1 + n * random_scales);
  diff_target[2] = mean(rowSums(1 - kappa)) - target_mean;
  i=2;
  while(T) {
    i = i+1;
    if(i > max_iter) {i = i-1; break;}
    log_scale[i] = log_scale[i-1] - diff_target[i-1]*(log_scale[i-1]-log_scale[i-2])/(diff_target[i-1]-diff_target[i-2]);
    random_scales = 1 / (slab_precision + 1/(exp(2*log_scale[i]) * lambda_sq));
    kappa = 1 / (1 + n * random_scales);
    diff_target[i] = mean(rowSums(1 - kappa)) - target_mean;
    if(abs(diff_target[i] - diff_target[i-1]) < tol) {break;}
  }

  return(list(scale = exp(log_scale[i]),
              diff_from_target = abs(diff_target[i]),
              iter = i,
              prior_num = diff_target[i] + target_mean));
}

