#' Matching Estimators for Multiple Treatments from Yang et al. (2016).
#'
#'
#' This function carries out matching on covariates or on propensity scores, and
#' is similar to the 'legacy' functions \code{\link{multilevelMatchX}} and
#' \code{\link{multilevelGPSMatch}}.
#'
#' @param Y A response vector (1 x n)
#' @param W A treatment vector (1 x n) with numerical values indicating
#'   treatment groups
#' @param X A covariate matrix (p x n) with no intercept. When
#'   match_on="existing", then X must be a vector (1 x n) of user-specified
#'   propensity scores.
#' @param M_matches Number of matches per unit for imputing potential outcomes,
#'   as in Abadie and Imbens (2006).
#' @param J_var_matches Number of matches when estimating \eqn{\sigma^2(X,W)} as
#'   in Abadie and Imbens (2006).
#' @param trimming an indicator of whether trimming the sample to ensure overlap
#' @param match_on User specifies "covariates" to match on raw covariates, or
#'   "existing" to match on user-supplied propensity score values, or "polr" or
#'   "multinom" to fit a propensity score model.
#' @param model_options A list of the options to pass to propensity model.
#'   Currently under development. Can only pass reference level to multinomial
#'   logistic regression.
#'
#' @return A list of output from \code{estimateTau}, including at most: \itemize{
#'
#'   \item \code{tauestimate}:  a vector of estimates for pairwise treatment
#'   effects
#'
#'   \item \code{varestimate}:  a vector of variance estimates for tauestimate,
#'   using Abadie and Imbens (2006)'s method
#'
#'   \item \code{varestimateAI2016}:  a vector of variance estimates for
#'   tauestimate, when matching on the generalized propensity score, using
#'   Abadie & Imbens (2016)'s method. This variance estimate takes into account
#'   of the uncertainty in estimating the GPS.
#'
#'   \item \code{analysis_idx}: a list containing the indices_kept (analyzed)
#'   and indices_dropped (trimmed) based on Crump et al. (2009)'s method.
#'
#'   }
#'
#' @seealso \code{\link{multilevelMatchX}}; \code{\link{multilevelGPSMatch}}
#'
#' @references Yang, S., Imbens G. W., Cui, Z., Faries, D. E., & Kadziola, Z.
#'   (2016) Propensity Score Matching and Subclassification in Observational
#'   Studies with Multi-Level Treatments. Biometrics, 72, 1055-1065.
#'   \url{https://doi.org/10.1111/biom.12505}
#'
#'   Abadie, A., & Imbens, G. W. (2006). Large sample properties of matching
#'   estimators for average treatment effects. econometrica, 74(1), 235-267.
#'   \url{https://doi.org/10.1111/j.1468-0262.2006.00655.x}
#'
#'   Abadie, A., & Imbens, G. W. (2016). Matching on the estimated propensity
#'   score. Econometrica, 84(2), 781-807.
#'   \url{https://doi.org/10.3982/ECTA11293}
#'
#'   Crump, R. K., Hotz, V. J., Imbens, G. W., & Mitnik, O. A. (2009). Dealing
#'   with limited overlap in estimation of average treatment effects.
#'   Biometrika, 96(1), 187-199. \url{https://doi.org/10.1093/biomet/asn055}
#'
#' @examples
#'   sim_data <- multilevelMatching::simulated_data
#'   Y <- sim_data$outcome
#'   W <- sim_data$treatment
#'   X <- as.matrix(sim_data[ ,-(1:2)])
#'   names(Y) <- paste0("ID", 1:length(Y))
#'
#'   trimming <- FALSE
#'   method <- c("covariates", "polr", "multinom")[2]
#'
#'   multiMatch(Y,W,X,trimming=trimming,match_on=method)
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

  ###                                          ###
  ### Step 1: Prepare data and other arguments ###
  ###                                          ###

  ## Argument checks and defensive programming
  prepared_data <- prepareData(
    Y = Y, W = W, X = X,
    match_on = match_on,
    trimming = trimming,
    model_options = model_options,
    M_matches = M_matches,
    J_var_matches = J_var_matches
  )
  X_covars <- prepared_data$X ## Necessary for AI2016 variance


  ## Fit treatment model to estimate propensity scores
  trt_model <- do.call(estimateTrtModel, prepared_data)
  prop_score_ests <- trt_model$prop_score_ests

  prepared_data <- append(
    prepared_data,
    list(
      M_matches = M_matches,
      J_var_matches = J_var_matches
    )
  )

  if (match_on != "covariates"){
    ## Instead of covariates, using propensity scores
    prepared_data$X <- prop_score_ests
  }


  ###                                          ###
  ### Step 2: Execute all matching procedures  ###
  ###                                          ###


  ## Carry out all estimation for all treatment levels
  matching_estimates <- do.call(matchAllTreatments, prepared_data)


  ###                                          ###
  ### Step 3: Organize into tidy, clean output ###
  ###                                          ###

  ## Tidy the results a little
  estimate_contrasts_args <- append(prepared_data, matching_estimates)
  results_list <- do.call(estimateTau, estimate_contrasts_args)
  tau_dfm <- results_list$tau_dfm

  impute_mat_sorted <- matching_estimates$Yiw
  impute_mat <- impute_mat_sorted[prepared_data$sorted_to_orig,]


  if (match_on == "multinom") {

  ###                                          ###
  ### Step 3B: Calculate AI2016 variance ests  ###
  ###                                          ###

    ## Collect a lot of objects from other function outputs
    est_var_AI2016_args <- append(estimate_contrasts_args, trt_model)
    est_var_AI2016_args$tau_dfm <- tau_dfm
    est_var_AI2016_args$X_covars <- X_covars

    ## Calculate AI2016 variance estimates for multinomial logistic regression
    est_var_AI2016 <- do.call(calcSigSqAI2016, est_var_AI2016_args)

    ## Add the AI2016 variance estimates to the tidy results dataframe
    tau_dfm$VarianceAI2016 <- est_var_AI2016$est_sigsq_AI2016
  }


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
    ## Additional information from the AI2016 variance estimation
    tidy_output$AI2016_args <- est_var_AI2016$AI2016_args
  }

  class(tidy_output) <- "multiMatch"
  tidy_output
}


# #' Calculate Variance as in Abadie & Imbens (2016) JASA
# #'
# #' This function calculates the estimated variance of matching estimator when
# #' matching on propensity scores that were estimated with multinomial logistic
# #' regression. This method takes into account the uncertainty from the treatment
# #' model. This function does no "heavy lifting" in that it only assembles pieces
# #' that have been estimated elsewhere (i.e., see \code{\link{estSigSq}} to see
# #' how the matching procedures are carried out).
# #'
# #' @inheritParams multiMatch
# #' @inheritParams matchAllTreatments
# #' @param X_covars The matrix of covariates (not the matrix of estimated GPS)
# #' @param prop_score_ests The matrix of propensity scores
# #' @param vcov_coeff The variance-covariance (inverse fisher information) matrix
# #'   from the multinomial logistic regression model
# #' @param tau_dfm The dataframe of estimates and estimated standard errors
# #'   output from \code{estimateTau}
# #' @param match_mat_AI2016 Information on matches from
# #'   \code{\link{matchAllTreatments}}
# #'
# #'
# #' @references Abadie, A., & Imbens, G. W. (2016). Matching on the estimated
# #'   propensity score. Econometrica, 84(2), 781-807.
# #'   \url{https://doi.org/10.3982/ECTA11293}
# #'
# #'
# #' @return A list with the updated \code{tau_dfm} including a column for
# #'   \code{VarianceAI2016}, and a list object \code{AI2016_args} with
# #'   potentially helpful extra information.
calcSigSqAI2016 <- function(
  W, X_covars, Y, N,
  num_trts, trt_levels,
  match_mat_AI2016,
  prop_score_ests,
  vcov_coeff,
  tau_dfm,
  ...
){
  est_sigsq_AI2016 <- rep(NA, NROW(tau_dfm))

  ncol_X <- dim(X_covars)[2]
  c_matrix <- matrix(0, nrow=N, ncol=(dim(X_covars)[2]+1)*(num_trts-1))
  c_vector <- matrix(0, nrow=num_trts, ncol=(ncol_X+1)*(num_trts-1))

  for ( trt_kk in 1:num_trts ) {
    col_name <- nameCols(trt_levels[trt_kk])
    Y11 <- matrix( Y[match_mat_AI2016[,c(col_name)]], ncol=2, byrow=FALSE)
    mY11 <- apply(Y11,1,mean)
    for ( trt_k in 1:(num_trts-1) ) {
      for ( trt_jj in 1:(ncol_X+1) ) {
        if ( trt_jj==1 ) {}
        if ( trt_jj>1 ) {
          X11 <- matrix( X_covars[match_mat_AI2016[,c(col_name)],(trt_jj-1)],
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

      est_sigsq_AI2016[result_row_num]  <-
        tau_dfm$Variance[result_row_num] -
        (
          t(c_vector[trt_jj,] + c_vector[trt_k,]) %*%
            vcov_coeff %*%
            ( c_vector[trt_jj,] + c_vector[trt_k,])
        )
    }
  }


  AI2016_args <- list(
    c_vector = c_vector,
    vcov_coeff = vcov_coeff,
    c_matrix = c_matrix
  )

  list(est_sigsq_AI2016 = est_sigsq_AI2016, AI2016_args = AI2016_args)
}
