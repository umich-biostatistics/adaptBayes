#' Numerical-based solution to the scale parameter c
#'
#' This function calculates a numerical-based solution to the scale parameter c in
#' the the equation three lines from the top of page 7 in Section 2 of Boonstra and
#' Barbaro. If desired, the user may request regional scale values for a partition
#' of the covariates into two regions, defined by the first npar1 covariates and the
#' second npar2 covariates, but this functionality was not used in Boonstra and Barbaro.
#'
#' @param target_mean1 (pos. reals): the desired prior number of effective parameters
#' (tilde xi_eff in Boonstra and Barbaro). If one scale parameter is desired, leave
#' target_mean2 = NA. An error will be thrown if target_mean1 > npar1 or if
#' target_mean2 > npar2.
#' @param target_mean2 (pos. reals): the desired prior number of effective parameters
#' (tilde xi_eff in Boonstra and Barbaro). If one scale parameter is desired, leave
#' target_mean2 = NA. An error will be thrown if target_mean1 > npar1 or if
#' target_mean2 > npar2.
#' @param npar1 (pos. integers): the number of covariates. If one scale parameter
#' is required, then leave npar2 = 0.
#' @param npar2 (pos. integers): the number of covariates. If one scale parameter
#' is required, then leave npar2 = 0.
#' @param local_dof (pos. integer) numbers indicating the degrees of freedom for
#' lambda_j and tau, respectively. Boonstra and Barbaro never considered
#' local_dof != 1 or global_dof != 1.
#' @param regional_dof (pos. integer) Not used in Boonstra and Barbaro.
#' @param global_dof (pos. integer) numbers indicating the degrees of freedom for
#' lambda_j and tau, respectively. Boonstra and Barbaro never considered
#' local_dof != 1 or global_dof != 1.
#' @param slab_precision (pos. real) the slab-part of the regularized horseshoe,
#' this is equivalent to (1/d)^2 in the notation of Boonstra and Barbaro
#' @param n (pos. integer) sample size
#' @param sigma (pos. real) square root of the assumed dispersion. In Boonstra and
#' Barbaro, this was always 2, corresponding to the maximum possible value:
#' sqrt(1/[0.5 * (1 - 0.5)]).
#' @param tol (pos. real) numerical tolerance for convergence of solution
#' @param max_iter (pos. integer) maximum number of iterations to run without
#' convergence before giving up
#' @param n_sim (pos. integer) number of simulated draws from the underlying student-t
#' hyperpriors to calculate the Monte Carlo-based approximation of the expectation.
#'
#'
#' @return A \code{list} containing the following named elements:
#' \itemize{
#'   \item{scale1}{}
#'   \item{diff_from_target1}{}
#'   \item{iter1}{}
#'   \item{prior_num1}{}
#'   \item{scale2}{}
#'   \item{diff_from_target2}{}
#'   \item{iter2}{}
#'   \item{prior_num2}{}
#' }
#'
#'
#' @export

solve_for_hiershrink_scale = function(target_mean1,
                                      target_mean2 = NA,
                                      npar1,
                                      npar2 = 0,
                                      local_dof = 1,
                                      regional_dof = -Inf,
                                      global_dof = 1,
                                      slab_precision = (1/15)^2,
                                      n,
                                      sigma = 2,
                                      tol = .Machine$double.eps^0.5,
                                      max_iter = 100,
                                      n_sim = 2e5
) {

  npar = npar1 + npar2;
  stopifnot(isTRUE(all.equal(npar%%1,0)));#Ensure integers
  do_local = (local_dof > 0);
  do_regional = (regional_dof > 0);
  do_global = (global_dof > 0);
  if(do_local) {
    lambda = matrix(rt(n_sim * npar,df = local_dof), nrow = n_sim);
  } else {
    lambda = matrix(1, nrow = n_sim, ncol = npar);
  }
  if(do_regional) {
    #Now do local-regional
    stopifnot(npar2 > 0 && isTRUE(all.equal(npar1%%1,0)) && isTRUE(all.equal(npar2%%1,0)));
    gamma_p = 1 + abs(rt(n_sim, df = regional_dof));
    gamma_q = 1 + abs(rt(n_sim, df = regional_dof));
    lambda[,1:npar1] = lambda[,1:npar1,drop=F] * gamma_p;
    lambda[,(npar1+1):(npar1+npar2)] = lambda[,(npar1+1):(npar1+npar2),drop=F] * gamma_q;
    rm(gamma_p,gamma_q);
  }
  if(do_global) {
    tau = rt(n_sim,df=global_dof);
    lambda[,1:npar1] = lambda[,1:npar1,drop=F] * tau;
    if(npar2 > 0) {
      lambda[,(npar1+1):(npar1+npar2)] = lambda[,(npar1+1):(npar1+npar2),drop=F] * tau;
    }
    rm(tau);
  }
  if(is.na(target_mean2)) {
    npar1 = npar;
  }
  stopifnot(target_mean1 > 0 && target_mean1 < npar1);#Ensure proper bounds
  log_scale1 = diff_target = numeric(max_iter);
  log_scale1[1] = log(target_mean1/(npar1 - target_mean1)*sigma/sqrt(n));
  random_scales = 1 / (slab_precision + 1/(exp(2*log_scale1[1]) * lambda^2));
  kappa = 1/(1+n*random_scales[,1:npar1,drop=F]/sigma^2);
  diff_target[1] = mean(rowSums(1-kappa)) - target_mean1;
  log_scale1[2] = 0.02 + log_scale1[1];
  random_scales = 1 / (slab_precision + 1/(exp(2*log_scale1[2]) * lambda^2));
  kappa = 1/(1+n*random_scales[,1:npar1,drop=F]/sigma^2);
  diff_target[2] = mean(rowSums(1-kappa)) - target_mean1;
  i=2;
  while(T) {
    i = i+1;
    if(i > max_iter) {i = i-1; break;}
    log_scale1[i] = log_scale1[i-1] - diff_target[i-1]*(log_scale1[i-1]-log_scale1[i-2])/(diff_target[i-1]-diff_target[i-2]);
    random_scales = 1 / (slab_precision + 1/(exp(2*log_scale1[i]) * lambda^2));
    kappa = 1/(1+n*random_scales[,1:npar1,drop=F]/sigma^2);
    diff_target[i] = mean(rowSums(1-kappa)) - target_mean1;
    if(abs(diff_target[i]-diff_target[i-1]) < tol) {break;}
  }
  scale1 = exp(log_scale1[i]);
  diff1 = abs(diff_target[i]);
  iter1 = i;
  prior_num1 = diff_target[i] + target_mean1;

  if(!is.na(target_mean2)) {
    stopifnot(target_mean2 > 0 && target_mean2 < npar2);
    log_scale2 = diff_target = numeric(max_iter);
    log_scale2[1] = log(target_mean2/(npar2 - target_mean2)*sigma/sqrt(n));
    random_scales = 1 / (slab_precision + 1/(exp(2*log_scale2[1]) * lambda^2));
    kappa = 1/(1+n*random_scales[,(npar1+1):(npar1+npar2),drop=F]/sigma^2);
    diff_target[1] = mean(rowSums(1-kappa)) - target_mean2;
    log_scale2[2] = 0.02 + log_scale2[1];
    random_scales = 1 / (slab_precision + 1/(exp(2*log_scale2[2]) * lambda^2));
    kappa = 1/(1+n*random_scales[,(npar1+1):(npar1+npar2),drop=F]/sigma^2);
    diff_target[2] = mean(rowSums(1-kappa)) - target_mean2;
    i=2;
    while(T) {
      i = i+1;
      if(i > max_iter) {i = i-1; break;}
      log_scale2[i] = log_scale2[i-1] - diff_target[i-1]*(log_scale2[i-1]-log_scale2[i-2])/(diff_target[i-1]-diff_target[i-2]);
      random_scales = 1 / (slab_precision + 1/(exp(2*log_scale2[i]) * lambda^2));
      kappa = 1/(1+n*random_scales[,(npar1+1):(npar1+npar2),drop=F]/sigma^2);
      diff_target[i] = mean(rowSums(1-kappa)) - target_mean2;
      if(abs(diff_target[i]-diff_target[i-1]) < tol) {break;}
    }
    scale2 = exp(log_scale2[i]);
    diff2 = abs(diff_target[i]);
    iter2 = i;
    prior_num2 = diff_target[i] + target_mean2;
  } else {
    scale2 = NA;
    diff2 = NA;
    iter2 = NA;
    prior_num2 = NA;
  }

  return(list(scale1 = scale1,
              diff_from_target1 = diff1,
              iter1 = iter1,
              prior_num1 = prior_num1,
              scale2 = scale2,
              diff_from_target2 = diff2,
              iter2 = iter2,
              prior_num2 = prior_num2));
}

