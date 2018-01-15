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
  model_options = list(reference_level = "1"),
  M_matches=1,
  J_var_matches=1
){

  warning("match_on='existing' needs further unit testing")

  ## Argument checks and defensive programming
  prepared_data <- prepareData(
    Y = Y, W = W, X = X,
    match_on = match_on,
    trimming = trimming,
    model_options = model_options
  )
  X <- prepared_data$X ## Necessary for AI2012 variance


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


  ## TODO: cleanup Abadie & Imbens 2012 variance estimation
  ##    migrate it to a subfunction
  W <- prepared_data$W
  # X <- prepared_data$X ## Do not re-assign, as it's been reassigned to PS
  Y <- prepared_data$Y
  N <- prepared_data$N
  num_trts <- prepared_data$num_trts
  trt_levels <- prepared_data$trt_levels
  N_per_trt <- prepared_data$N_per_trt
  num_contrasts <- prepared_data$num_contrasts
  analysis_idx <- prepared_data$analysis_idx

  if (match_on == "multinom") {
    match_mat_AI2012 <- matching_estimates$match_mat_AI2012

    tau_dfm$VarianceAI2012 <- NA

    vcov_coeff <- trt_model$vcov_coeff
    # I_inverse <- trt_model$vcov_coeff
    ## Adjustment term c'(I^-1)c
    X <- as.matrix(X)
    c_matrix <- matrix(0,N,(dim(X)[2]+1)*(num_trts-1))
    c_vector <- matrix(0,num_trts,(dim(X)[2]+1)*(num_trts-1))

    for(trt_kk in 1:num_trts){
      col_name <- nameCols(trt_levels[trt_kk])
      Y11 <- matrix( Y[match_mat_AI2012[,c(col_name)]], ncol=2, byrow=FALSE)
      mY11 <- apply(Y11,1,mean)
      for(trt_k in 1:(num_trts-1)){
        for(trt_jj in 1:(dim(X)[2]+1)){
          if(trt_jj==1){}
          if(trt_jj>1){
            X11 <- matrix( X[match_mat_AI2012[,c(col_name)],(trt_jj-1)],
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
            c_matrix[,(dim(X)[2]+1)*(trt_k-1)+trt_jj] <- c1_X1Y
          }
        }
      }
      c_vector[trt_kk,] <- apply(c_matrix,2,mean)
    }

    for(trt_jj in 1:(num_trts-1)){
      for(trt_k in (trt_jj+1):num_trts){
        result_row_num <-
          which( tau_dfm$Param == nameContrast(
            trt1 = trt_levels[trt_jj], trt2 = trt_levels[trt_k]))

        tau_dfm$VarianceAI2012[result_row_num]  <-
          tau_dfm$Variance[result_row_num] -
          t(c_vector[trt_jj,] + c_vector[trt_k,]) %*% vcov_coeff %*% ( c_vector[trt_jj,] + c_vector[trt_k,])
      }
    }
    AI2012_args <- list(
      c_vector = c_vector,
      vcov_coeff = vcov_coeff,
      c_matrix = c_matrix
    )
  }



  tidy_output <- list(
    results = tau_dfm,
    analysis_idx = analysis_idx,
    mu = results_list$mu_dfm,
    impute_mat = matching_estimates$Yiw[prepared_data$sorted_to_orig,],
    estimate_args = estimate_contrasts_args,
    model = trt_model$prop_score_model,
    propensity_scores = prop_score_ests,
    impute_match_data = matching_estimates$impute_match_data
  )

  if (match_on == "multinom") {
    tidy_output$AI2012_args <- AI2012_args
  }

  tidy_output
}
