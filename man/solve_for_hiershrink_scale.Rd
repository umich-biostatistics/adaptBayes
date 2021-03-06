% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solve_for_hiershrink_scale.R
\name{solve_for_hiershrink_scale}
\alias{solve_for_hiershrink_scale}
\title{Numerical-based solution to the scale parameter c}
\usage{
solve_for_hiershrink_scale(
  target_mean1,
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
  n_sim = 2e+05
)
}
\arguments{
\item{target_mean1}{(pos. reals): the desired prior number of effective parameters
(tilde xi_eff in Boonstra and Barbaro). If one scale parameter is desired, leave
target_mean2 = NA. An error will be thrown if target_mean1 > npar1 or if
target_mean2 > npar2.}

\item{target_mean2}{(pos. reals): the desired prior number of effective parameters
(tilde xi_eff in Boonstra and Barbaro). If one scale parameter is desired, leave
target_mean2 = NA. An error will be thrown if target_mean1 > npar1 or if
target_mean2 > npar2.}

\item{npar1}{(pos. integers): the number of covariates. If one scale parameter
is required, then leave npar2 = 0.}

\item{npar2}{(pos. integers): the number of covariates. If one scale parameter
is required, then leave npar2 = 0.}

\item{local_dof}{(pos. integer) numbers indicating the degrees of freedom for
lambda_j and tau, respectively. Boonstra and Barbaro never considered
local_dof != 1 or global_dof != 1.}

\item{regional_dof}{(pos. integer) Not used in Boonstra and Barbaro.}

\item{global_dof}{(pos. integer) numbers indicating the degrees of freedom for
lambda_j and tau, respectively. Boonstra and Barbaro never considered
local_dof != 1 or global_dof != 1.}

\item{slab_precision}{(pos. real) the slab-part of the regularized horseshoe,
this is equivalent to (1/d)^2 in the notation of Boonstra and Barbaro}

\item{n}{(pos. integer) sample size}

\item{sigma}{(pos. real) square root of the assumed dispersion. In Boonstra and
Barbaro, this was always 2, corresponding to the maximum possible value:
sqrt(1/[0.5 * (1 - 0.5)]).}

\item{tol}{(pos. real) numerical tolerance for convergence of solution}

\item{max_iter}{(pos. integer) maximum number of iterations to run without
convergence before giving up}

\item{n_sim}{(pos. integer) number of simulated draws from the underlying student-t
hyperpriors to calculate the Monte Carlo-based approximation of the expectation.}
}
\value{
A \code{list} containing the following named elements:
\itemize{
  \item{scale1}{}
  \item{diff_from_target1}{}
  \item{iter1}{}
  \item{prior_num1}{}
  \item{scale2}{}
  \item{diff_from_target2}{}
  \item{iter2}{}
  \item{prior_num2}{}
}
}
\description{
This function calculates a numerical-based solution to the scale parameter c in
the the equation three lines from the top of page 7 in Section 2 of Boonstra and
Barbaro. If desired, the user may request regional scale values for a partition
of the covariates into two regions, defined by the first npar1 covariates and the
second npar2 covariates, but this functionality was not used in Boonstra and Barbaro.
}
