
#' Perform matching methods to estimate (or impute) potential outcomes for one
#' treatment level (level kk)
#'
#' This function executes the matching methods to estimate the components for
#' the point estimates of the population-level estimands. The components will be
#' assembled into the estimator in the function \code{\link{estimateTau}}.
#'
#' @inheritParams matchAllTreatments
#' @param est_options Options for carrying out matching for point estimation.
#' @param vec_same_trt Boolean vector indicating which individuals observed to
#'   have treatment level kk
#' @param which_same_trt Numeric vector providing indices of individuals
#'   observed to have treatment level kk
#' @param vec_diff_trt Boolean vector indicating which individuals observed to
#'   have treatment level NOT equal to kk
#' @param which_diff_trt Numeric vector providing indices of individuals
#'   observed to have treatment level  NOT equal to kk
#'
#' @return A list, many of which are necessary arguments for
#'   \code{\link{estimateTau}}: \itemize{ \item \code{Yiw} The fully-imputed
#'   data set of all potential outcomes. \item \code{mean_Yiw} Mean of Yiw
#'   across all units i. These are the estimates of the population-level mean
#'   estimands. See \code{\link{estimateTau}}. \item \code{Kiw} The vector of
#'   number of times unit i used as a match in imputation. \item
#'   \code{impute_match_data} Extra information from the main matching
#'   procedure. }
matchImputePO <- function(
  X, W, Y,
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

  #                                    #
  # Step 3: Calculate some components  #
  #                                    #

  diff_trt_match_data <- diff_trt_match$mdata

  K_factor <- factor(
    x = diff_trt_match$index.control,
    levels = which_same_trt
  )
  Kiw_kk <- table(K_factor)

  Yiw_kk <- rep(NA, length(Y))
  ## TODO: allow for M>1 number of matches
  Yiw_kk[which_same_trt] <- Y[which_same_trt]
  Yiw_kk[which_diff_trt] <-
    diff_trt_match_data$Y[which(diff_trt_match_data$Tr==0)]

  mean_po_kk  <- stats::weighted.mean(
    x = c(Yiw_kk[which_same_trt], Yiw_kk[which_diff_trt]),
    w = c(rep(1,length(which_same_trt)), diff_trt_match$weights)
  )


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

