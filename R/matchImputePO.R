
# #' Perform matching methods to estimate (or impute) potential outcomes for one
# #' treatment level (level \code{kk})
# #'
# #' This function executes the matching methods to estimate the components for
# #' the point estimates of the population-level estimands. The components will be
# #' assembled into the estimator in the function \code{\link{estimateTau}}.
# #'
# #' @inheritParams matchAllTreatments
# #' @param est_options Options for carrying out matching for point estimation.
# #' @param vec_same_trt Boolean vector indicating which individuals observed to
# #'   have treatment level kk
# #' @param which_same_trt Numeric vector providing indices of individuals
# #'   observed to have treatment level kk
# #' @param vec_diff_trt Boolean vector indicating which individuals observed to
# #'   have treatment level NOT equal to kk
# #' @param which_diff_trt Numeric vector providing indices of individuals
# #'   observed to have treatment level  NOT equal to kk
# #'
# #' @return A list, many of which are necessary arguments for
# #'   \code{\link{estimateTau}}: \itemize{ \item \code{Yiw} The fully-imputed
# #'   data set of all potential outcomes. \item \code{mean_Yiw} Mean of Yiw
# #'   across all units i. These are the estimates of the population-level mean
# #'   estimands. See \code{\link{estimateTau}}. \item \code{Kiw} The vector of
# #'   number of times unit i used as a match in imputation. \item
# #'   \code{impute_match_data} Extra information from the main matching
# #'   procedure. }
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
matchImputePO <- function(
  X,
  # W,
  Y,
  vec_diff_trt, which_diff_trt,
  vec_same_trt, which_same_trt,
  est_options
){

  #                                    #
  # Step 1: Prepare matching arguments #
  #                                    #

  diff_trt_matching_args <- append(
    est_options,
    list(
      Y = Y,
      Tr = vec_diff_trt,
      X = X
    )
  )

  #                                    #
  # Step 2: Call matching procedure    #
  #                                    #

  diff_trt_match <- do.call(Matching::Match, diff_trt_matching_args)
  ## Defensive programming / CYA
  if ( any(diff_trt_match$weights ==0) ||
       !is.null(diff_trt_match$index.dropped)  ){
    warn_message <- paste(
      "Matching::Match() has dropped units;",
      "this may result in errors or unintentional consequences in",
      "multilevelMatching::multiMatch()."
    )
    warning(warn_message)
  }


  #                                    #
  # Step 3: Calculate some components  #
  #                                    #

  est_kk_list <- wrangleImputations(
    match_output = diff_trt_match,
    M_matches = diff_trt_matching_args$M,
    Y = Y,
    which_same_trt = which_same_trt,
    which_diff_trt = which_diff_trt
  )
}


# #' Calculate imputed potential outcomes and other estimation components
# #'
# #' This is a plumbing function, called from \code{\link{matchImputePO}}; they
# #' output the same list of information.
# #'
# #' @inheritParams multiMatch
# #' @inheritParams matchImputePO
# #' @param match_output The output from \code{\link[Matching]{Match}} when
# #'   matching different treatment levels to impute potential outcomes (for point
# #'   estimation of causal estimands).
# #'
# #' @return A list, many of which are necessary arguments for
# #'   \code{\link{estimateTau}}: \itemize{ \item \code{Yiw} The fully-imputed
# #'   data set of all potential outcomes. \item \code{mean_Yiw} Mean of
# #'   \code{Yiw} across all units i. These are the estimates of the
# #'   population-level mean estimands. See \code{\link{estimateTau}}. \item
# #'   \code{Kiw} The vector of number of times unit i used as a match in
# #'   imputation. \item \code{impute_match_data} Extra information from the main
# #'   matching procedure. }
# #'
wrangleImputations <- function(
  match_output,
  M_matches,
  Y,
  which_same_trt,
  which_diff_trt
){

  diff_trt_match_data <- match_output$mdata


  K_factor <- factor(
    x = match_output$index.control,
    levels = which_same_trt
  )
  Kiw_kk <- table(K_factor)

  Yiw_kk <- rep(NA, length(Y))

  orig_outcomes <- Y[which_same_trt]
  ## Wrangle outcomes in case there is one-to-many matching (M>1)
  outcomes_list <- averageMultipleMatches(
    num_matches = M_matches,
    matched_outcomes =
      diff_trt_match_data$Y[which(diff_trt_match_data$Tr==0)]
  )
  matched_outcomes <- outcomes_list$matched_outcomes

  Yiw_kk[which_same_trt] <- orig_outcomes
  Yiw_kk[which_diff_trt] <- matched_outcomes

  mean_po_kk <- mean(Yiw_kk)

  ## Save below for later; this may be helpful if matching drops units and some
  ##  weights are equal to 0
  #
  #mean_po_kk  <- stats::weighted.mean( x =
  #c(Yiw_kk[which_same_trt],
  #diff_trt_match_data$Y[which(diff_trt_match_data$Tr==0)]), w =
  #c(rep(1,length(which_same_trt)), diff_trt_match$weights) )

  stopifnot(!any(is.na(Yiw_kk)))
  #                                    #
  # Step 4: Organize Output            #
  #                                    #

  list_out <- list(
    mean_po_kk = mean_po_kk,
    Kiw_kk = Kiw_kk,
    Yiw_kk = Yiw_kk,
    impute_match_data_kk = diff_trt_match_data
  )
}

# #' Plumbing function for one-to-many matches
# #'
# #' This is called from \code{\link{wrangleImputations}} and from
# #' \code{\link{calcSigSqAI2006}}.
# #'
# #' @param num_matches Either \code{M_matches} or \code{J_var_matches}
# #' @param matched_outcomes These are the \code{num_matches}-many imputed
# #'   potential outcomes. These need to be averaged if \code{num_matches}>1.
# #' @param orig_outcomes When called from \code{\link{calcSigSqAI2006}}, these
# #'   are repeated outcomes that simply need to be subsetted; in the case of
# #'   \code{\link{wrangleImputations}} this is left to \code{NULL} and is
# #'   ignored.
# #'
# #' @return A list including the averaged \code{matched_outcomes}, and also
# #'   \code{orig_outcomes} which is sometimes \code{NULL}.
# #'
averageMultipleMatches <- function(
  num_matches, matched_outcomes, orig_outcomes = NULL
){

  if (num_matches==1){
    orig_outcomes <- orig_outcomes
    matched_outcomes <- matched_outcomes
  } else {
    if (!is.null(orig_outcomes)){
      unique_indivs <- seq(1, length(orig_outcomes), by = num_matches)
      orig_outcomes <- orig_outcomes[unique_indivs]
    }
    matched_outcomes_matrix <- matrix(
      matched_outcomes,
      ncol = (1/num_matches) * length(matched_outcomes),
      nrow = num_matches,
      byrow = FALSE
    )
    matched_outcomes <- (1/num_matches) * colSums(matched_outcomes_matrix)
  }


  list(
    orig_outcomes = orig_outcomes,
    matched_outcomes = matched_outcomes
  )
}
