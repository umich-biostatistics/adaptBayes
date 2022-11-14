#' Projection matrix (or list of projection matrices)
#'
#'
#' Helper function to make projection matrix (or list of projection matrices), which
#' is P in the notation of Boonstra and Barbaro.
#'
#'
#'
#' @param x_curr_orig (matrices) matrices with equal numbers of rows and p & q columns,
#' respectively. These are used to estimate the joint association between the original
#' and augmented covariates, which is needed for the imputation model
#' @param x_curr_aug (matrices) matrices with equal numbers of rows and p & q columns,
#' respectively. These are used to estimate the joint association between the original
#' and augmented covariates, which is needed for the imputation model
#' @param eigenvec_hist_var (matrix) pxp matrix with each row corresponding to an
#' eigenvector. This is v_0 in Boonstra and Barbaro.
#' @param imputes_list (list) list of length-2 vectors, with each vector containing
#' the lower and upper indices of the imputations to use for a projection matrix
#' in the SAB method. This is best explained with an example: if
#' imputes_list = list(c(1,15),c(16,100),c(1,100)), then three projection matrices
#' will be returned. One will be based upon the first 15 imputations from a run of
#' MICE, the second based upon the last 85 imputations from that same run (i.e. the 16th to100th
#' imputations), and the third will be based upon all 100 imputations from this same run.
#' This is coded as such to allow for flexible exploration of the impact of number of
#' imputations or variability due to imputations.
#' @param seed_start (pos. integer) random seed to start each imputation
#' @param predictorMatrix (matrix) (p+q)x(p+q) matrix equivalent to argument of the
#' same name in in the 'mice()' function (type '?mice'). It is specially calculated
#' based upon a monotone missingness pattern (x^o is fully observed, x^a is not) Thus
#' it samples from
#'
#'
#'
#' @return A list as long as the the length of imputes_list, with each element containing
#' a different projection matrix using the indices of the imputations specified in the
#' corresponding element of imputes_list.
#'
#' @importFrom stats cor
#' @export

create_projection = function(x_curr_orig,
                             x_curr_aug,
                             eigenvec_hist_var,
                             imputes_list = list(c(1,500)),
                             seed_start = sample(.Machine$integer.max,1),
                             predictorMatrix = NULL) {

  p = ncol(x_curr_orig);
  q = ncol(x_curr_aug);
  orig_covariates = colnames(x_curr_orig);
  aug_covariates = colnames(x_curr_aug);
  x_all = cbind(x_curr_orig,x_curr_aug);
  stopifnot(class(imputes_list) == "list");

  n_imputes = max(unlist(lapply(imputes_list,max)));

  if(is.null(predictorMatrix)) {
    predictorMatrix11 = matrix(0,nrow = p,ncol = p);
    predictorMatrix12 = matrix(0,nrow = p,ncol = q);
    predictorMatrix21 = matrix(1,nrow = q,ncol = p);
    predictorMatrix22 = matrix(1,nrow = q,ncol = q);
    predictorMatrix22[upper.tri(predictorMatrix22,diag = T)] = 0;
    predictorMatrix = rbind(cbind(predictorMatrix11,predictorMatrix12),
                            cbind(predictorMatrix21,predictorMatrix22));
    colnames(predictorMatrix) = rownames(predictorMatrix) = c(orig_covariates,aug_covariates);
    rm(predictorMatrix11,predictorMatrix12,predictorMatrix21,predictorMatrix22);
  } else {
    stopifnot(dim(predictorMatrix) == c(p+q,p+q));
  }

  #Create data to be marginalized
  #Column 1 identifies each eigenvector, including an additional row for the intercept offset;
  #Columns 2:(p+1) are original covariates;
  #Columns (p+2):(p+q+1) are the augmented covariates to be imputed and averaged;
  dat_for_marg = cbind(c(1:p,nrow(x_all)+1),
                       rbind(eigenvec_hist_var,0),
                       matrix(NA,p+1,q));
  colnames(dat_for_marg) = c("ID",orig_covariates,aug_covariates);
  dat_for_marg = data.matrix(dat_for_marg);
  #Store one copy for each of the different imputations to use
  dat_for_marg = lapply(1:length(imputes_list), function(x) dat_for_marg);

  #Now impute the augmented covariates based upon the artificially created original covariates, i.e. the eigenvectors, using
  #the empiric associations in current data.
  curr_impute = mice::mice(rbind(x_all,
                                 #original covariates are constant across the different imputations, so we can just use the first
                                 dat_for_marg[[1]][,c(orig_covariates,aug_covariates)]),
                           printFlag = F,
                           predictorMatrix = predictorMatrix,#user-provided, or if missing, created above
                           m = n_imputes,#
                           maxit = 1,#only one iteration is needed, due to the monotone missingness pattern
                           method = "pmm",
                           seed = seed_start);
  curr_impute = data.matrix(mice::complete(curr_impute,"long"));
  for(k in 1:(p+1)) {
    for(j in 1:length(imputes_list)) {
      #Fill in the data_for marg
      dat_for_marg[[j]][dat_for_marg[[j]][,"ID"] == dat_for_marg[[j]][k,"ID"], aug_covariates] =
        colMeans(curr_impute[which(curr_impute[,".id"] == k + nrow(x_all))[imputes_list[[j]][1]:imputes_list[[j]][2]],aug_covariates,drop=F]);
    }
  }

  #If there is any perfect correlation between predictors, the imputer will return NA. In this case, we fill in the missing data
  #with whatever variable was perfectly correlated with it
  for(j in 1:length(imputes_list)) {
    while(any(is.na(dat_for_marg[[j]][,aug_covariates]))) {
      correlated_column = which(colSums(is.na(dat_for_marg[[j]][,aug_covariates]))>0)[1];
      dat_for_marg[[j]][,names(correlated_column)] =  dat_for_marg[[j]][,setdiff(colnames(x_all)[abs(cor(x_all,x_curr_aug[,correlated_column]) - 1) < sqrt(.Machine$double.eps)],names(correlated_column))[1]];
    }
  }
  #v0_inv is constant across the different imputation sets, so we can just use the first set
  v0_inv = solve(dat_for_marg[[1]][1:p,orig_covariates,drop=F]);
  projections = vector("list",length(imputes_list));
  for(j in 1:length(imputes_list)) {
    projections[[j]] = v0_inv %*% (dat_for_marg[[j]][1:p,aug_covariates,drop=F] - dat_for_marg[[j]][rep(p+1,p),aug_covariates,drop=F]);
  }
  projections;
}
