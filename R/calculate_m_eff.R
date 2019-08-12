
#DESCRIPTION: This function is the inverse of 'solve_for_hiershrink_scale'. Instead of providing a desired effective number of 
#parameters, the user provides the scale value(s), which is c in the notation of Boonstra and Barbaro, and the the function gives
#the implied prior number of effective parameters based upon this. As with 'solve_for_hiershrink_scale', the user can provide
#one global scale parameter (scale1, leaving scale2 = NA) that applies to all parameters, or two regional scale parameters (scale1, 
#scale2), that applies to a partition of the parameters as defined by the first npar1 parameters and the second npar2 parameters.
#
#
#ARGUMENTS:
#target_mean1, target_mean2 (pos. reals): the desired prior number of effective parameters (tilde xi_eff in Boonstra and Barbaro). 
#If one scale parameter is desired, leave target_mean2 = NA. An error will be thrown if target_mean1 > npar1 or if 
#target_mean2 > npar2. 

#' Compute inverse of 'solve_for_hiershrink_scale'
#'
#' Instead of providing a desired effective number of parameters, the user provides 
#' the scale value(s), which is c in the notation of Boonstra and Barbaro, and the the 
#' function gives the implied prior number of effective parameters based upon this. 
#' As with 'solve_for_hiershrink_scale', the user can provide one global scale parameter 
#' (scale1, leaving scale2 = NA) that applies to all parameters, or two regional scale 
#' parameters (scale1, scale2), that applies to a partition of the parameters as 
#' defined by the first npar1 parameters and the second npar2 parameters.
#'
#' @param scale1 
#' @param scale2 
#' @param npar1 
#' @param npar2 
#' @param local_dof 
#' @param regional_dof 
#' @param global_dof 
#' @param slab_precision 
#' @param n 
#' @param sigma 
#' @param tol 
#' @param max_iter 
#' @param n_sim 
#'
#' @return list containing prior numbers 1 and 2
#'

calculate_m_eff = function(scale1,
                           scale2 = NA,
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
  do_local = (local_dof > 0);
  do_regional =(regional_dof > 0);
  do_global = (global_dof > 0);
  npar = npar1 + npar2;
  stopifnot(isTRUE(all.equal(npar%%1,0)));#Ensure integers
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
  random_scales = 1 / (slab_precision + 1/(scale1^2*lambda^2));
  if(is.na(scale2)) {
    npar1 = npar;
  }
  kappa1 = 1/(1+n*random_scales[,1:npar1,drop=F]/sigma^2);
  prior_num1 = mean(rowSums(1-kappa1))
  
  if(!is.na(scale2)) {
    random_scales = 1 / (slab_precision + 1/(scale2^2*lambda^2));
    kappa2 = 1/(1+n*random_scales[,(npar1+1):(npar1+npar2),drop=F]/sigma^2);
    prior_num2  = mean(rowSums(1-kappa2));
  } else {
    prior_num2 = NA;
  }
  
  return(list(prior_num1 = prior_num1,
              prior_num2 = prior_num2));
}
