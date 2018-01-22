#' Matching Estimators for Mutltiple Yreatments from Yang et al (2016).
#'
#'
#' This function carries out matching on covariates or on propensity scores, and
#' is similar to the 'legacy' functions \code{\link{multilevelMatchX}} and
#' \code{\link{multilevelGPSMatch}}.
#'
#' @param Y a continuous response vector (1 x n)
#' @param W a treatment vector (1 x n) with numerical values indicating
#'   treatment groups
#' @param X A covariate matrix (p x n) with no intercept. When
#'   match_on="existing", then X must be a vector (1 x n) of user-specified
#'   propensity scores.
#' @param M_matches Number of matches per unit for imputing potential outcomes,
#'   as in Abadie and Imbens 2006. Currently can only support M=1.
#' @param J_var_matches Number of matches when estimating sigmasq(X,W) as in
#'   Abadie and Imbens 2006. Currently can only support J=1.
#' @param trimming an indicator of whether trimming the sample to ensure overlap
#' @param match_on "multinom", "polr", "existing", or "covariates",
#' @param model_options A list of the options to pass to propensity model.
#'   Currently under development. Can only pass reference level to multinomial
#'   logisitc regression.
#'
#' @return according to \code{\link{estimateTau}}, including at most: \itemize{
#'
#'   \item tauestimate:  a vector of estimates for pairwise treatment effects
#'
#'   \item varestimate:  a vector of variance estimates for tauestimate, using
#'   Abadie&Imbens(2006)'s method
#'
#'   \item varestimateAI2012:  a vector of variance estimates for tauestimate,
#'   when matching on the generalized propensity score, using
#'   Abadie&Imbens(2012)'s method. This variance estimate takes into account of
#'   the uncertainty in estimating the GPS.
#'
#'   \item analysis_idx: a list containing the indices_kept (analyzed) and
#'   indices_dropped (trimmed) based on Crump et al. (2009)'s method.
#'
#'   }
#'
#' @seealso \code{\link{multilevelMatchX}}; \code{\link{multilevelGPSMatch}}
#'
#' @examples
#'   X <- as.matrix(c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9),ncol=1)
#'   Y <- c(102,105,120,130,100,80,94,108,96)
#'   W <- c(1,1,1,3,2,3,2,1,2)
#'   multiMatch(Y,W,X,trimming=0,match_on="multinom")
#'   multiMatch(Y,W,X,trimming=1,match_on="multinom")
#'   multiMatch(Y,W,X,trimming=0,match_on="polr")
#'   multiMatch(Y,W,X,trimming=0,match_on="covariates")
#'
#'
#' @export
multiMatch <- function(
  Y,W,X,
  trimming = NULL,
  match_on,
  model_options = list(reference_level = sort(W)[1]),
  M_matches=1,
  J_var_matches=1
){


  ## Argument checks and defensive programming
  prepared_data <- prepareData(
    Y = Y, W = W, X = X,
    match_on = match_on,
    trimming = trimming,
    model_options = model_options
  )
  X_covars <- prepared_data$X ## Necessary for AI2012 variance


  ## Fit treatment model to estimate propensity scores
  trt_model <- do.call(estimateTrtModel, prepared_data)
  prop_score_ests <- trt_model$prop_score_ests

  prepared_data <- append(
    prepared_data,
    list(
      M_matches = M_matches,
      J_var_matches = J_var_matches
    ))

  if (match_on != "covariates"){
    ## Instead of covariates, using propensity scores
    prepared_data$X <- prop_score_ests
  }



  ## Carry out all estimation for all treatment levels
  matching_estimates <- do.call(matchAllTreatments,prepared_data)

  ## Tidy the results a little
  estimate_contrasts_args <- append(prepared_data, matching_estimates)
  results_list <- do.call(estimateTau,estimate_contrasts_args)
  tau_dfm <- results_list$tau_dfm


  if (match_on == "multinom") {
    est_var_AI2012_args <- append(estimate_contrasts_args,trt_model)
    est_var_AI2012_args$tau_dfm <- tau_dfm
    est_var_AI2012_args$X_covars <- X_covars

    ## Estimate AI2012 variance for multinomial logistic regression
    est_var_AI2012 <- do.call(estVarAI2012, est_var_AI2012_args)

    ##This will be added to the output. So will AI2012_args.
    tau_dfm <- est_var_AI2012$tau_dfm
  }

  impute_mat_sorted <- matching_estimates$Yiw
  impute_mat <- impute_mat_sorted[prepared_data$sorted_to_orig,]

  tidy_output <- list(
    results = tau_dfm,
    analysis_idx = prepared_data$analysis_idx,
    mu = results_list$mu_dfm,
    impute_mat = impute_mat,
    estimate_args = estimate_contrasts_args,
    model = trt_model$prop_score_model,
    propensity_scores = prop_score_ests,
    impute_match_data = matching_estimates$impute_match_data,
    impute_mat_sorted = impute_mat_sorted
  )

  if (match_on == "multinom") {
    tidy_output$AI2012_args <- est_var_AI2012$AI2012_args
  }

  tidy_output
}


#' Estimate Variance as in Abadie & Imebens 2012 JASA
#'
#' This function estimates variance of matching estimator when matching on
#' propensity scores that were estimated with multinomial logistic regression.
#' This method takes into account the uncertainty from the treatment model.
#'
#' @inheritParams multiMatch
#' @inheritParams matchAllTreatments
#' @param X_covars The matrix of covariates (not the matrix of estimated GPS)
#' @param prop_score_ests The matrix of propensity scores
#' @param vcov_coeff The variance-covariance (inverse fisher information) matrix
#'   from the multinomial logistic regression model
#' @param tau_dfm The dataframe of estimates and estimated standard errors
#'   output from \code{\link{estimateTau}}
#' @param match_mat_AI2012 Information on matches from
#'   \code{\link{matchAllTreatments}}
#'
#' @return A list with the updated \code{tau_dfm} includeing a column for
#'   \code{VarianceAI2012}, and a list object \code{AI2012_args} with
#'   potentially helpful extra information.
estVarAI2012 <- function(
  W, X_covars, Y, N,
  num_trts, trt_levels,
  match_mat_AI2012,
  prop_score_ests,
  vcov_coeff,
  tau_dfm,
  ...
){
  tau_dfm$VarianceAI2012 <- NA

  ncol_X <- dim(X_covars)[2]
  c_matrix <- matrix(0, nrow=N, ncol=(dim(X_covars)[2]+1)*(num_trts-1))
  c_vector <- matrix(0, nrow=num_trts, ncol=(ncol_X+1)*(num_trts-1))

  for ( trt_kk in 1:num_trts ) {
    col_name <- nameCols(trt_levels[trt_kk])
    Y11 <- matrix( Y[match_mat_AI2012[,c(col_name)]], ncol=2, byrow=FALSE)
    mY11 <- apply(Y11,1,mean)
    for ( trt_k in 1:(num_trts-1) ) {
      for ( trt_jj in 1:(ncol_X+1) ) {
        if ( trt_jj==1 ) {}
        if ( trt_jj>1 ) {
          X11 <- matrix( X_covars[match_mat_AI2012[,c(col_name)],(trt_jj-1)],
                         ncol=2, byrow=FALSE)
          mX11 <- apply(X11,1,mean)

          c1_X1Y_temp <- apply((X11-mX11)*(Y11-mY11),1,sum)
          if ( trt_kk == (trt_k+1) ) {
            mult_factor <- (1-prop_score_ests[,trt_k+1])
          } else
            if ( trt_kk != (trt_k+1) ) {
              mult_factor <- (-1)*prop_score_ests[,trt_k+1]
            }
          c1_X1Y <- c1_X1Y_temp*mult_factor
          c_matrix[,(ncol_X+1)*(trt_k-1)+trt_jj] <- c1_X1Y
        }
      }
    }
    c_vector[trt_kk,] <- apply(c_matrix,2,mean)
  }

  for ( trt_jj in 1:(num_trts-1) ) {
    for ( trt_k in (trt_jj+1):num_trts ) {
      result_row_num <-
        which( tau_dfm$Param == nameContrast(
          trt1 = trt_levels[trt_jj], trt2 = trt_levels[trt_k]))

      tau_dfm$VarianceAI2012[result_row_num]  <-
        tau_dfm$Variance[result_row_num] -
        (
          t(c_vector[trt_jj,] + c_vector[trt_k,]) %*%
            vcov_coeff %*%
            ( c_vector[trt_jj,] + c_vector[trt_k,])
        )
    }
  }


  AI2012_args <- list(
    c_vector = c_vector,
    vcov_coeff = vcov_coeff,
    c_matrix = c_matrix
  )

  list(tau_dfm = tau_dfm, AI2012_args = AI2012_args)
}
