
#' Estimate Imputed Outcomes and AI06 Variance Components
#'
#' Mostly calls a subfunction for each treatment level, which calls subfunctions
#' that either carry out the matched imputation or the variance estimation.
#'
#' @inheritParams multiMatch
#' @inheritParams estimateTrtModel
#' @param num_trts The number of unique treatments (3+), a scalar
#' @param trt_levels A vector (of length \code{num_trts} providing the unique
#'   treatment levels
#' @param N_per_trt A vector (of length \code{num_trts}) indicating the number
#'   of units observed to have each treatment level
#' @param N The total number of units
#' @param unit_ids_sorted The unit identifiers in correct, sorted order.
#'
#' @return A list including at most: \itemize{
#'
#'   \item Yiw: A matrix of all imputed (or observed) potential outcomes for
#'   each unit \item mean_Yiw: A vector of the average (across all units) of
#'   estimated/imputed potential outcomes \item sigsqiw: Estimated variance
#'   components (from AI2006) for each unit \item Kiw: The vector of number of
#'   times unit i used as a match \item impute_match_data: extra information
#'   from the main matching procedure \item match_mat_AI2012: When
#'   \code{match_on=`multinom`} this additional information will be output for
#'   \code{\link{estVarAI2012}} }
#'
#'   A few of the necessary arguments are output from the
#'   \code{\link{reorderByTreatment}} function.
#'
matchAllTreatments <- function(
  N, X, W, Y,
  num_trts, trt_levels, N_per_trt, unit_ids_sorted,
  M_matches,  J_var_matches,
  match_on,
  ...
){

  blank_mat <- matrix(NA, ncol=num_trts, nrow=N)
  colnames(blank_mat) <- trt_levels
  rownames(blank_mat) <- unit_ids_sorted

  list_out <- list(
    impute_match_data = list(), ## extra information from the main matching proc
    sigsqiw = blank_mat[,1, drop=FALSE], # The estimated sigma squared for each unit i,
    Kiw = blank_mat[,1, drop=FALSE],   # The vector of number of times unit i used as a match,
    Yiw = blank_mat, # The full imputed data set,
    mean_Yiw = blank_mat[1,]) # Mean of Yiw across all units i

  basic_matching_args <- list(
    distance.tolerance = 0,
    ties = FALSE, ## Ties will be randomly broken.
    Weight = 2) ## Mahalanobis

  diff_trt_matching_basic_args <- append( ## For imputing potential outcomes
    basic_matching_args,
    list(X = X, Y = Y, M = M_matches))

  same_trt_match_basic_args <- append( ## For estimating AI2006 variance
    basic_matching_args,
    list(M = J_var_matches))

  match_once_basic_args <- list( ## For carrying out both matching procedures
    X = X, W = W, Y = Y,
    diff_trt_matching_basic_args = diff_trt_matching_basic_args,
    same_trt_match_basic_args = same_trt_match_basic_args,
    match_on = match_on)
    # GPSM = GPSM)


  if (match_on == "multinom") {
  # if (GPSM=="multinomiallogisticReg") {

    match_mat_AI2012 <- matrix(NA,N,num_trts*2)
    colnames(match_mat_AI2012) <- nameCols(trt_levels)
    list_out$match_mat_AI2012 <- match_mat_AI2012

    inside_match_AI2012_basic_args <- list(
      M = 1,#M_matches,
      distance.tolerance=0,
      ties=FALSE,
      Weight=2
    )
    match_once_basic_args$inside_match_AI2012_basic_args <-
      inside_match_AI2012_basic_args
    match_once_basic_args$match_mat_AI2012 <- match_mat_AI2012
  }

  for(kk in 1:num_trts){

    this_trt_level <- trt_levels[kk]
    N_this_trt  <- N_per_trt[kk]
    indices_to_skip <- ifelse( kk==1, 0, sum(N_per_trt[1:(kk-1)]) )
    trtd_indiv_indices <- (1:N_per_trt[kk]) + indices_to_skip

    if (match_on %in% c("multinom", "polr", "existing")) {
      X_GPS_kk <- X[,kk]
      match_once_basic_args$diff_trt_matching_basic_args$X <- X_GPS_kk
      stopifnot(!is.matrix(match_once_basic_args$diff_trt_matching_basic_args$X))
      match_once_basic_args$X <- X_GPS_kk
    }
    match_once_args <- append(
      match_once_basic_args,
      list(
        this_trt_level = this_trt_level,
        trtd_indiv_indices = trtd_indiv_indices,
        N_this_trt = N_this_trt))

    if (match_on == "multinom") {
      # if (GPSM=="multinomiallogisticReg"),
      match_once_args$match_mat_col_name <- nameCols(this_trt_level)
    }

    ## Carry out imputation-matching and variance estimation for kk^th treatment
    match_trt_kk <- do.call(matchOneTreatment, match_once_args)

    list_out$sigsqiw[which(W==trt_levels[kk]), 1] <- match_trt_kk$sigsqiw_kk
    list_out$Kiw[which(W==trt_levels[kk]), 1] <-  match_trt_kk$Kiw_kk
    list_out$mean_Yiw[kk] <-  match_trt_kk$mean_po_kk
    list_out$Yiw[,kk] <- match_trt_kk$Yiw_kk
    list_out$impute_match_data[[kk]] <- match_trt_kk$impute_match_data_kk
    if (match_on == "multinom"){
      list_out$match_mat_AI2012[,2*kk+(-1:0)] <-
        match_trt_kk$match_mat_AI2012_kk[,2*kk+(-1:0)]
    }
  }

  list_out
}



# #' performs imputation of potential outcomes, and estimates AI06 standard
# #' errors, for one treatment level
matchOneTreatment <- function(
  X, W, Y,
  trtd_indiv_indices, this_trt_level, N_this_trt,
  diff_trt_matching_basic_args, same_trt_match_basic_args,
  ## For AI2012 matching on multinomial logistic reg
  inside_match_AI2012_basic_args = NULL,
  match_mat_AI2012=NULL, match_mat_col_name = NULL,
  match_on
){

  ## Helpers for determining treatment levels
  vec_diff_trt <- W != this_trt_level
  vec_same_trt <- W == this_trt_level
  which_diff_trt <- which(vec_diff_trt)
  which_same_trt <- which(vec_same_trt)

  shared_matching_args <- list(X = X, Y = Y, which_same_trt = which_same_trt)

  match_imputation_args <- append(
    shared_matching_args,
    list(
      vec_diff_trt = vec_diff_trt,
      vec_same_trt = vec_same_trt,
      which_diff_trt = which_diff_trt,
      trtd_indiv_indices = trtd_indiv_indices,
      diff_trt_matching_basic_args = diff_trt_matching_basic_args))

  match_sigsqiw_args <- append(
    shared_matching_args,
    list(
      N_this_trt = N_this_trt,
      match_on = match_on,
      ## For AI2012 MLR
      match_mat_AI2012 = match_mat_AI2012,
      inside_match_AI2012_basic_args = inside_match_AI2012_basic_args,
      match_mat_col_name = match_mat_col_name,
      vec_diff_trt = vec_diff_trt,
      ## For all
      same_trt_match_basic_args = same_trt_match_basic_args))



  ## First, diff-matching for imputing potential outcomes
  imputed_pos_kk <- do.call(matchImputePO, match_imputation_args)
  ## TODO:: check if colmean(Yiw)==mean_Yiw


  ## Second, matching within-treatment to get the AI06 (& AI2012) variance estimate
  sigsqiw_kk_est <- do.call(estSigSq,match_sigsqiw_args)

  ## Add variance to estimate output
  # imputed_pos_kk$sigsqiw_kk <- sigsqiw_kk_est
  imputed_pos_kk <- append(
    imputed_pos_kk,
    sigsqiw_kk_est
  )

  imputed_pos_kk
}








# ## carries out matching for imputing potential outcomes from one treatment to
# ## individuals on all other treatment levels
matchImputePO <- function(
  X, W, Y,
  vec_diff_trt, which_diff_trt,
  vec_same_trt, which_same_trt,
  trtd_indiv_indices,
  diff_trt_matching_basic_args
){

  diff_trt_matching_kk_args <- list(Tr = vec_diff_trt)
  diff_trt_matching_args <- append(
    diff_trt_matching_basic_args,diff_trt_matching_kk_args
  )
  diff_trt_match <- do.call(Matching::Match,diff_trt_matching_args)
  diff_trt_match_data <- diff_trt_match$mdata

  weighted_mean_args <- list(
    x = c(Y[which_same_trt], diff_trt_match_data$Y[which(diff_trt_match_data$Tr==0)]),
    w = c(rep(1,length(which_same_trt)), diff_trt_match$weights)
  )

  ## Matching to impute potential outcomes
  mean_po_kk  <- do.call(stats::weighted.mean, weighted_mean_args)

  Kiw_kk <- table(factor(diff_trt_match$index.control,levels=trtd_indiv_indices))
  Yiw_kk <- rep(NA, length(Y))
  ## TODO: allow for M>1 number of matches
  Yiw_kk[which_same_trt] <- Y[which_same_trt]
  Yiw_kk[which_diff_trt] <- diff_trt_match_data$Y[which(diff_trt_match_data$Tr==0)]

  list_out <- list(
    mean_po_kk = mean_po_kk,
    Kiw_kk = Kiw_kk,
    Yiw_kk = Yiw_kk,
    impute_match_data_kk = diff_trt_match_data
  )
}


# #' Estimates sigma_squared by the matching procedure as in Abadie&Imbens2006
estSigSq <- function(
  X, Y,
  which_same_trt, N_this_trt,
  same_trt_match_basic_args,
  match_on,
  ## these are for AI2012
  vec_diff_trt,
  inside_match_AI2012_basic_args,
  match_mat_col_name,
  match_mat_AI2012

){

  if ( !is.matrix(X) ) {
    X_mat_same_trt <- as.matrix(X[which_same_trt])
  } else {
    X_mat_same_trt <- as.matrix(X[which_same_trt, ])
  }

  ## repeat trtd individuals 2x
  outcome_repeated <- rep(Y[which_same_trt], times=2)
  trt_repeated <- rep(c(1,0), each=N_this_trt)
  restriction_matrix <- matrix(
    c(1:(2*N_this_trt),
      rep(-1,N_this_trt)),
    nrow = N_this_trt, ncol = 3, byrow = FALSE
  ) ## this restricts so we can't match an individual to itself

  same_trt_match_kk_args <- list(
    Y  = outcome_repeated,
    Tr = trt_repeated,
    X  = rbind(X_mat_same_trt, X_mat_same_trt),
    restrict = restriction_matrix
  )

  same_trt_match_args <- append(
    same_trt_match_basic_args,
    same_trt_match_kk_args
  )

  ## Matching to estimate variance component as in Abadie&Imbens2006
  same_trt_match <- do.call(Matching::Match, args = same_trt_match_args)

  same_trt_match_data <- same_trt_match$mdata
  J_var_matches <- same_trt_match_basic_args$M
  ## TODO: take mean of the J matches when J>1
  sigsq_est_temp <- ((J_var_matches)/(1+J_var_matches))*
    (same_trt_match_data$Y[which(same_trt_match_data$Tr==1)] -
       same_trt_match_data$Y[which(same_trt_match_data$Tr==0)])^2


  out_sigma <- list(
    sigsqiw_kk = sigsq_est_temp
  )

  ## Abadie&Imbens2012 variance estimator for multinomial logistic regression
  # if (GPSM=="multinomiallogisticReg") {
  if (match_on=="multinom") {


    # col_name <- nameCols(thistrt) ## match_mat_col_name


    # find two outsiders closest
    outside_match_args <- list(
      Y = Y,
      Tr = vec_diff_trt,
      X = X,
      distance.tolerance = 0,
      ties = FALSE,
      Weight = 2,
      M = 2
    )
    outside_match <- do.call(Matching::Match, outside_match_args)
    # findmatch1 <- Matching::Match(Y=Y,Tr=vec_diff_trt,X=PF.fit[,kk],distance.tolerance=0,ties=FALSE,Weight=2,M=2)
    match_mat_AI2012[unique(outside_match$index.treated),match_mat_col_name] <-
      matrix(outside_match$index.control,ncol=2,byrow=TRUE)


    X_GPS_vec_same_trt <- X[which_same_trt]

    # find one insider closest
    inside_match_kk_args <- list(
      Y = outcome_repeated, #rep(Y[which_same_trt],times=2),
      Tr = rep(c(0,1), each=N_this_trt),
      X = c(X_GPS_vec_same_trt,X_GPS_vec_same_trt),
      restrict = restriction_matrix
    )
    inside_match_args <- append(
      inside_match_AI2012_basic_args,inside_match_kk_args
    )

    inside_match <- do.call(Matching::Match, inside_match_args)

    match_mat_AI2012[which_same_trt,match_mat_col_name] <-
      matrix( c(which_same_trt, which_same_trt[inside_match$index.control]),
              ncol=2, byrow=FALSE)

    out_sigma$match_mat_AI2012_kk <- match_mat_AI2012
  }

  out_sigma
}
