

# #' Perform matching methods to estimate variance for one treatment level
# #'
# #' This function executes the matching methods to estimate the components for
# #' estimating the variance (sigma squared) of the matching estimator. All
# #' matching methods (specified with \code{match_on}) will estimate the variance
# #' as in Abadie and Imbens (2006). When \code{match_on = "multinom"} this function
# #' will carry out additional matching procedures to estimate the variance as
# #' described in Abadie and Imbens (2016).
# #'
# #' @inheritParams multiMatch
# #' @inheritParams matchAllTreatments
# #' @inheritParams matchImputePO
# #' @param var_options Options for carrying out matching for variance estimation.
# #' @param N_this_trt The number of units observed to have this treatment level.
# #' @param var_args_AI2016 A list of arguments for carrying out matching
# #'   procedures to estimate components in the \code{VarianceAI2016} variance estimates
# #'   (see \code{\link{calcSigSqAI2016}}).
# #'
# #'   Note that these variance components are combined in
# #'   \code{\link{estimateTau}} (and perhaps \code{\link{calcSigSqAI2016}}).
# #'
# #' @return A list of one or two elements. \code{sigsqiw_kk} is a vector with the
# #'   estimated variance component for units observed to have the \code{kk}^th treatment
# #'   level. When \code{match_on = "multinom"}, the list will also have an
# #'   element for \code{match_mat_AI2016_kk_two_cols}, which is a \code{Ntot}-by-2
# #'   matrix with matching information that will eventually be passed to
# #'   \code{\link{calcSigSqAI2016}}.
# #'
# #' @references Abadie, A., & Imbens, G. W. (2006). Large sample properties of
# #'   matching estimators for average treatment effects. econometrica, 74(1),
# #'   235-267. \url{https://doi.org/10.1111/j.1468-0262.2006.00655.x}
# #'
# #'   Abadie, A., & Imbens, G. W. (2016). Matching on the estimated propensity
# #'   score. Econometrica, 84(2), 781-807.
# #'   \url{https://doi.org/10.3982/ECTA11293}
estSigSq <- function(
  X, Y,
  which_same_trt, N_this_trt,
  var_options,
  match_on,
  var_args_AI2016
){

  #                                    #
  # Step 1: Prepare matching arguments #
  #                                    #

  if ( !is.matrix(X) ) {
    X_mat_same_trt <- as.matrix(X[which_same_trt])
  } else {
    X_mat_same_trt <- as.matrix(X[which_same_trt, ])
  }

  ## Repeat trtd individuals 2x
  outcome_repeated <- rep(Y[which_same_trt], times=2)
  trt_repeated <- rep(c(1,0), each=N_this_trt)
  restriction_matrix <- matrix(
    c(1:(2*N_this_trt),
      rep(-1,N_this_trt)),
    nrow = N_this_trt, ncol = 3, byrow = FALSE
  ) ## restriction_matrix will not allow an individual to match to itself

  same_trt_match_args <- append(
    var_options,
    list(
      Y  = outcome_repeated,
      Tr = trt_repeated,
      X  = rbind(X_mat_same_trt, X_mat_same_trt),
      restrict = restriction_matrix
    )
  )

  #                                    #
  # Step 2: Call matching procedure    #
  #                                    #

  ## Matching to estimate variance component as in Abadie&Imbens2006
  same_trt_match <- do.call(Matching::Match, args = same_trt_match_args)


  #                                    #
  # Step 3: Calculate some components  #
  #                                    #

  out_sigma <- list(
    sigsqiw_kk = calcSigSqAI2006(
      match_output = same_trt_match,
      J_var_matches = same_trt_match_args$M
      ## originally J_var_matches in `multiMatch()`
    )
  )

  ## Abadie & Imbens (2016) variance estimator for multinomial logistic regression
  if (match_on=="multinom") {

    ## First, matching to find two outsiders (different treatments) closest

    outside_match_args <- append(
      var_args_AI2016$var_options_AI2016,
      list(
        Y = Y,
        Tr = var_args_AI2016$vec_diff_trt, ## Outsiders (between-trt-levels)
        X = X,
        M = 2 ## Two matches per unit
      )
    )

    outside_match <- do.call(Matching::Match, outside_match_args)


    ## Then, match to find one insider (same treatment) closest

    X_GPS_vec_same_trt <- X[which_same_trt]
    inside_match_args <- append(
      var_args_AI2016$var_options_AI2016,
      list(
        Y = outcome_repeated,
        Tr = rep(c(0,1), each=N_this_trt),
        X = c(X_GPS_vec_same_trt,X_GPS_vec_same_trt), ## Insiders (same treatment)
        restrict = restriction_matrix, ## Don't match to self
        M = 1 ## Only one
      )
    )

    inside_match <- do.call(Matching::Match, inside_match_args)

    #                                    #
    # Step 4b: Organize output, again    #
    #                                    #

    mat_kk <-
      var_args_AI2016$match_mat_AI2016_kk_two_cols

    mat_kk[unique(outside_match$index.treated), ] <-
      matrix( outside_match$index.control, ncol=2, byrow=TRUE)

    mat_kk[which_same_trt,] <-
      matrix( c(which_same_trt, which_same_trt[inside_match$index.control]),
              ncol=2, byrow=FALSE)



    ## Add match_mat to output of estSigSq()
    out_sigma$match_mat_AI2016_kk_two_cols <- mat_kk
  }

  out_sigma
}


# #' Calculates \code{sigsqiw}
# #'
# #' This function estimates the conditional variance as seen in equation 14 in
# #' Section 4.1 of Abadie and Imbens 2006 Econometrica. The matching procedure in
# #' \code{\link{estSigSq}} matches units within the same treatment level (level
# #' kk), and compares their outcomes to estimate a variance component. This
# #' function was introduced to extend this variance matching procedure to
# #' \code{J_var_matches >=1}, as it takes care of some of the bookkeeping aspects
# #' of one-to-many matching.
# #'
# #' @inheritParams multiMatch
# #' @param match_output Output of the \code{Matching::Match()} function for the
# #'   same-treatment matching, from \code{\link{estSigSq}}.
# #'
# #' @return A vector of the \code{sigsqiw} values for those individuals observed
# #'   to have the \code{kk}th treatment level
calcSigSqAI2006 <- function(match_output,J_var_matches){
  md <- match_output$mdata
  J_factor <- ( (J_var_matches)/(1+J_var_matches) )

  outcomes_list <- averageMultipleMatches(
    num_matches = J_var_matches,
    orig_outcomes = md$Y[which(md$Tr==1)],
    matched_outcomes = md$Y[which(md$Tr==0)]
  )
  orig_outcomes <- outcomes_list$orig_outcomes
  matched_outcomes <- outcomes_list$matched_outcomes

  ## The estimated conditional variance as in AI2006
  ## for those units observed to have the kk^th treatment level
  sigsqiw_kk <- J_factor * ( orig_outcomes - matched_outcomes )^2
}

