
# #' Estimate Imputed Outcomes and AI06 Variance Components
# #'
# #' Mostly calls a subfunction for each treatment level, which calls subfunctions
# #' that either carry out the matched imputation or the variance estimation.
# #'
# #' @inheritParams multiMatch
# #' @inheritParams estimateTrtModel
# #' @param num_trts The number of unique treatments (3+), a scalar
# #' @param trt_levels A vector (of length \code{num_trts} providing the unique
# #'   treatment levels
# #' @param N_per_trt A vector (of length \code{num_trts}) indicating the number
# #'   of units observed to have each treatment level
# #' @param N The total number of units
# #' @param unit_ids_sorted The unit identifiers in correct, sorted order.
# #'
# #' @return A list including at most: \itemize{
# #'
# #'   \item \code{Yiw}: A matrix of all imputed (or observed) potential outcomes
# #'   for each unit \item \code{mean_Yiw}: A vector of the average (across all
# #'   units) of estimated/imputed potential outcomes \item \code{sigsqiw}:
# #'   Estimated variance components (from Abadie and Imbens (2006)'s method) for
# #'   each unit \item \code{Kiw}: The vector of number of times unit i used as a
# #'   match \item \code{impute_match_data}: extra information from the main
# #'   matching procedure \item \code{match_mat_AI2016}: When
# #'   \code{match_on=`multinom`} this additional information will be output for
# #'   \code{\link{calcSigSqAI2016}} (from Abadie and Imbens (2016)'s method) }
# #'
# #'   A few of the necessary arguments are output from the
# #'   \code{\link{reorderByTreatment}} function.
# #'
# #'
# #' @references Yang, S., Imbens G. W., Cui, Z., Faries, D. E., & Kadziola, Z.
# #'   (2016) Propensity Score Matching and Subclassification in Observational
# #'   Studies with Multi-Level Treatments. Biometrics, 72, 1055-1065.
# #'   \url{https://doi.org/10.1111/biom.12505}
# #'
# #'   Abadie, A., & Imbens, G. W. (2006). Large sample properties of
# #'   matching estimators for average treatment effects. econometrica, 74(1),
# #'   235-267. \url{https://doi.org/10.1111/j.1468-0262.2006.00655.x}
# #'
# #'   Abadie, A., & Imbens, G. W. (2016). Matching on the estimated propensity
# #'   score. Econometrica, 84(2), 781-807.
# #'   \url{https://doi.org/10.3982/ECTA11293}
# #'
matchAllTreatments <- function(
  N, X, W, Y,
  num_trts, trt_levels, N_per_trt, unit_ids_sorted,
  M_matches,  J_var_matches,
  match_on,
  ...
){

  ##                                                            ##
  ## Step 1: Creating objects that will hold the important data ##
  ##                                                            ##


  blank_mat <- matrix(NA, ncol=num_trts, nrow=N)
  colnames(blank_mat) <- trt_levels
  rownames(blank_mat) <- unit_ids_sorted

  list_out <- list(
    # The fully-imputed data set of all potential outcomes,
    Yiw = blank_mat,
    # Mean of Yiw across all units i
    mean_Yiw = blank_mat[1,],
    # The vector of number of times unit i used as a match in imputation
    Kiw = blank_mat[,1, drop=FALSE],
    # Extra information from the main matching procedure
    impute_match_data = list(),
    # The estimated sigma squared for each unit i,
    sigsqiw = blank_mat[,1, drop=FALSE]
  )
  colnames(list_out$sigsqiw) <- "sigsqiw"
  colnames(list_out$Kiw) <- "Kiw"


  ##                                                            ##
  ## Step 2: Preparing basic arguments invariant to treatments  ##
  ##                                                            ##

  basic_matching_args <- list(
    distance.tolerance = 0,
    ties = FALSE, ## Ties will be randomly broken.
    Weight = 2 ## Mahalanobis
  )

  est_options <- append( ## For imputing potential outcomes
    basic_matching_args, list(M = M_matches) ## TODO: Allow for M>1.
  )

  var_options <- append( ## For estimating AI2006 variance
    basic_matching_args, list(M = J_var_matches)  ## TODO: Allow for J>1.
  )

  if (match_on == "multinom") {
    var_args_AI2016 <- list(var_options_AI2016 = basic_matching_args)

    match_mat_AI2016 <- matrix(NA,N,num_trts*2)
    colnames(match_mat_AI2016) <- nameCols(trt_levels)
    ## This matrix is the important output for AI2016 variance estimation
    list_out$match_mat_AI2016 <- match_mat_AI2016
  }

  ## Hacky but useful for incorporating matching on GPS (see below).
  if (match_on == "covariates") { X_kk <- X }


  ##                                                             ##
  ## Step 3: Loop over each treatment level and execute matching ##
  ##                                                             ##

  for(kk in 1:num_trts){

    ## Helpers for determining treatment levels, etc
    this_trt_level <- trt_levels[kk]
    N_this_trt  <- N_per_trt[kk]
    vec_diff_trt <- W != this_trt_level
    vec_same_trt <- W == this_trt_level
    which_diff_trt <- which(vec_diff_trt)
    which_same_trt <- which(vec_same_trt)

    if (match_on %in% c("multinom", "polr", "existing")) {
      ## Matching on only the kk^th column of propensity scores for trt kk.
      X_kk <- X[,kk, drop = TRUE]
      stopifnot(!is.matrix(X_kk))
    }

    ## These used in both matchImputePO() and estSigSq()
    shared_args <- list(
      X = X_kk, Y = Y, which_same_trt = which_same_trt
    )


    ##                                                             ##
    ## Step 4a: Matching for point estimation with matchImputePO() ##
    ##                                                             ##

    ##### matchImputePO()
    ##### First, matching between distinct treatment levels
    #####   for imputing potential outcomes

    est_args <- append(
      shared_args,
      list(
        vec_diff_trt = vec_diff_trt,
        vec_same_trt = vec_same_trt,
        which_diff_trt = which_diff_trt,
        est_options = est_options
      )
    )

    est_kk_list <- do.call(matchImputePO, est_args)
    ## TODO:: testthat colmean(Yiw)==mean_Yiw, or similar calculation

    ## Add mean ests stuff to output
    list_out$Kiw[which(W==trt_levels[kk]), 1] <-  est_kk_list$Kiw_kk
    list_out$mean_Yiw[kk] <-  est_kk_list$mean_po_kk
    list_out$Yiw[,kk] <- est_kk_list$Yiw_kk
    list_out$impute_match_data[[kk]] <- est_kk_list$impute_match_data_kk

    ##                                                           ##
    ## Step 4b: Matching for variance estimation with estSigSq() ##
    ##                                                           ##

    ##### estSigSq()
    ##### Second, matching within-treatment
    #####   to get the AI06 (& AI2016) variance estimate

    var_args <- append(
      shared_args,
      list(
        N_this_trt = N_this_trt,
        match_on = match_on,
        var_options = var_options
      )
    )

    if (match_on == "multinom") {
      var_args_AI2016$match_mat_AI2016_kk_two_cols <-
        match_mat_AI2016[,2*kk+(-1:0)]

      var_args_AI2016$vec_diff_trt <- vec_diff_trt

      ## Additional information for estimating AI2016 variance
      var_args$var_args_AI2016 <- var_args_AI2016
    }

    sigsqiw_kk_est <- do.call(estSigSq,var_args)


    ## Add variance ests stuff to output
    list_out$sigsqiw[which(W==trt_levels[kk]), 1] <- sigsqiw_kk_est$sigsqiw_kk
    if (match_on == "multinom"){
      list_out$match_mat_AI2016[,2*kk+(-1:0)] <-
        sigsqiw_kk_est$match_mat_AI2016_kk_two_cols
    }

  }
  ## End of the treatment level loop


  list_out
}



