

#' Prepare data for estimation
#'
#' A series of checks, tests, re-ordering, and other operations to prepare the
#' data for matching. This function can be run standalone, before running
#' \code{\link{multiMatch}}.
#'
#' @inheritParams multiMatch
#'
#'
#' @return A list of information, including the \code{X, W, Y} arguments after
#'   sorting observeations, and information on \code{unit_ids}, etc.
#'
#' @examples
#'
#'  sim_data <- multilevelMatching::simulated_data
#'  Y <- sim_data$outcome
#'  W <- sim_data$treatment
#'  X <- as.matrix(sim_data[ ,-(1:2)])
#'  names(Y) <- paste0("ID", 1:length(Y))
#'
#'  trimming <- FALSE
#'  method <- c("covariates", "polr", "multinom")[2]
#'
#'  prepared_data <- prepareData(
#'    Y = Y,
#'    W = W,
#'    X = X,
#'    match_on = "polr",
#'    trimming = FALSE,
#'    model_options = list(reference_level = sort(W)[1]),
#'    M_matches = 3,
#'    J_var_matches = 2
#'  )
#' @export
prepareData <- function(
  Y, W, X, match_on, trimming=NULL, model_options,
  M_matches, J_var_matches
){

  if ( (M_matches<=0)|| (M_matches %%1 !=0) ) {
    stop("M_matches must be a positive integer")
  }
  if ( (J_var_matches<=0)|| (J_var_matches %%1 !=0) ) {
    stop("J_var_matches must be a positive integer")
  }


  data_ids <- determineIDs(Y=Y, W=W, X=X)
  Y <- data_ids$Y
  W <- data_ids$W
  X <- data_ids$X
  unit_ids <- unit_ids_orig <- data_ids$unit_ids


  if (!is.null(model_options)){
    if (!is.list(model_options)){
      stop("model_options must be a list or NULL")
    }
    if (match_on == "multinom"){
      if (!"reference_level" %in% names(model_options)){
        stop("User must supply model_options$reference_level for multinomial model")
      } ## defensive programming for correct level variable?
    }
  }

  if ((length(class(X))!=1 && class(X)!="matrix") || class(X)!="matrix") {
    long_message <- "Dataframes, lists, and some other types are not supported at this time"
    warning("It is recommended that X is a matrix. X will be coerced to a matrix.")
    if (is.list(X) || is.data.frame(X)) {stop(long_message)}
    if (is.factor(X)) {
      warning("It is recommended that X is a matrix.\n Factor variables are not recommended at this time")
    }

    X <- try(as.matrix(X))
    if ("try-error" %in% class(X)){
      stop("Error with coercing X to matrix. Please specify X as a matrix")
    }

  }

  ## Defensive checks
  argChecks(Y=Y, W=W, X=X, N=NULL)


  if (!( match_on %in% c("covariates", "polr", "multinom", "existing") )) {
    stop('match_on must be "polr", "multinom", "existing" or "covariates"')
  }



  if (is.null(trimming)){
    if (match_on == "multinom"){
      stop("trimming is a necessary argument when match_on='multinom': set trimming to FALSE or TRUE")
    } else {
      trimming <- 0
    }
  }

  if (trimming==0) {
    analysis_idx <- NULL
    # analysis_idx <- list(ids_kept = unit_ids)
    } else
    if (trimming ==1){
      message_no_trimming <- "trimming is currently only supported for multinomial logistic regression"
      if (match_on == "covariates") {
        stop(message_no_trimming)
      } else if (match_on == "existing") {
        stop(message_no_trimming)
      } else if (match_on == "polr") {
        stop(message_no_trimming)
      } else if (match_on == "multinom"){

        #PF modeling
        W.ref <- stats::relevel(as.factor(W),ref=1)
        temp <- utils::capture.output(PF.out <- nnet::multinom(W.ref~X))
        ## TODO: make this into a subfunction?
        PF.fit <- stats::fitted(PF.out)

        ## identify sufficient overlap
        overlap_idx <- overlap(PF.fit)$idx
        W <- W[overlap_idx]
        X <- X[overlap_idx,,drop=FALSE]
        Y <- Y[overlap_idx]

        ## organize indices into 'kept' and 'dropped' for clean output
        all_indices <- 1:length(Y)
        dropped_idx <- all_indices[!all_indices%in% overlap_idx]
        unit_ids <- unit_ids_orig[overlap_idx]
        analysis_idx <- list(
          indices_kept = overlap_idx,
          indices_dropped = dropped_idx,
          ids_kept = unit_ids,
          ids_dropped = unit_ids_orig[dropped_idx]
        )
      }
    } else {stop("trimming must be set to FALSE or TRUE")}



  ## Order the treatment increasingly
  ordered_data_list <- reorderByTreatment(W=W,X=X,Y=Y, unit_ids_unsorted = unit_ids)

  if ( any(ordered_data_list$N_per_trt <= J_var_matches) ) {
    stop("J_var_matches is too large: there are not enough treatments in at least one treatment level to successfully match")
  }
  if ( any(ordered_data_list$N_per_trt < M_matches) ) {
    stop("M_matches is too large: there are not enough treatments in at least one treatment level to successfully match")
  }

  if ( match_on == "existing" ) {
    ## Check that "existing" X has the correctly specified number of levels
    #Warnings fixed with GH issue BarkleyBG/multilevelMatching/#4 ## see NEWS
    #for v0.2.4. warning("match_on='existing' needs further unit testing")
    #warning("user-supplied propensity scores has not been passing checks as of
    #2017-03-26")
    if (nrow(X) != ordered_data_list$N) {
      stop("user-supplied propensity scores (through argument X and
           match_on='existing') should be a matrix with number of rows
           equal to length of Y and W")
    }
    if (ncol(X) != ordered_data_list$num_trts) {
      stop("user-supplied propensity scores (through argument X and
           match_on='existing') should be a matrix with number of columns
           equal to the number of treatment levels in W")
    }

    ## TODO: Check the row sum of X is 1 - debug1
    ##
    indiv_ps_sums <- rowSums(X)
    if (any(indiv_ps_sums!=1)) {
      stop(
        paste0(
          "All rows of X must sum to 1. The following rows do not
          (output limited to 5 row numbers): c(",
          paste( which(indiv_ps_sums!=1)[1:5],collapse = ","),
          ")"
        )
      )
    }
  }

  ordered_data_list$analysis_idx <- analysis_idx
  ordered_data_list$match_on <- match_on
  ordered_data_list$model_options <- model_options

  ordered_data_list
}

# #' Determines Unique Unit Identifiers
# #'
# #' This function attempts to determine unique identifying information,
# #' \code{unit_ids}, for each unit in the dataset. Users can apply this function
# #' on their raw data ahead of using \code{\link{multiMatch}} to ensure that the
# #' matching procedure will work. \code{unit_ids} will be used to identify study
# #' units in some of the information output from \code{\link{multiMatch}}.
# #'
# #' @inheritParams multiMatch
# #'
# #' @return \code{unit_ids}
# #'
determineIDs <- function(Y,W,X){


  names_list <- list(
    "Y" = getIDs(Y),
    "W" = getIDs(W),
    "X" = getIDs(X)
  )

  has_null_names <- unlist(lapply(names_list, function(x){
    ifelse(is.null(x), 1, 0)}))

  num_units <- ifelse(is.matrix(W) || is.data.frame(W), NROW(W), length(W))

  num_null_names <- sum(has_null_names)
  which_not_null <- which(!has_null_names )


  has_names <- paste( c("Y", "W", "X")[which_not_null] ,
                      collapse = " and ")
  unit_id_stop_message <- paste(
    "Non-identical names were specified for arguments:",
    has_names,
    ". Please specify identical names/rownames to Y, W and X.",
    "These are considered to be the unit IDs."
  )

  if ( num_null_names == 3 ) {
    message("no unit IDs supplied; unit_ids will be assigned generically")
    unit_ids <- paste0("unit_", 1:num_units)
  } else
    if ( num_null_names == 2 ) {
      unit_ids <- names_list[[which_not_null]]
    } else
      if ( num_null_names == 1 ) {
        names1 <- names_list[[which_not_null[1]]]
        names2 <- names_list[[which_not_null[2]]]

        if (identical(names1, names2)) {
          unit_ids <- names1
        } else {
          stop(unit_id_stop_message)
        }
      } else
        if ( num_null_names == 0 ){
          if (
            identical(names_list[[1]], names_list[[2]]) &&
            identical(names_list[[1]], names_list[[3]]) &&
            identical(names_list[[3]], names_list[[2]])
          ) {
            unit_ids <- names_list[[1]]
          } else {
            stop(unit_id_stop_message)
          }

        } else {
          stop("Error with unit_ids. Please assign names/rownames of Y,W, and X to be the same")
        }

  if ( any(duplicated(unit_ids)) ) {
    stop("unit_ids are duplicated. Please assign unique names/rownames of Y,W and X")
  }

  list(
    Y = setIDs(Y,unit_ids),
    W = setIDs(W,unit_ids),
    X = setIDs(X,unit_ids),
    unit_ids = unit_ids
  )
}

# #' Get/grab identifiers from vector/matrix/dataframe
# #'
# #' This is a helper function for \code{\link{determineIDs}}
# #'
# #' @param x An object
# #'
# #' @return \code{names(x)}, \code{row.names(x)}, or \code{NULL}.
getIDs <- function(x){
  if (is.vector(x)) {names <- names(x)} else
    if (is.matrix(x)) {names <- rownames(x)} else
      if (is.data.frame(x)) {names <- row.names(x)} else {
        stop(
          paste0(
            "Y,W, and X must be vectors, matrices, or perhaps data.frames. ",
            "Factors should be avoided."
          )
        )
      }
  names
}

# #' Set identifiers from vector/matrix/dataframe
# #'
# #' This is a helper function for \code{\link{determineIDs}}
# #'
# #' @param x An object
# #' @param unit_ids Character vector to identify study units.
# #'
# #' @return The object \code{x} after setting its \code{names(x)},
# #'   \code{row.names(x)}, or \code{rownames(x)} to \code{unit_ids}.
setIDs <- function(x,unit_ids){
  if (is.vector(x)) {names(x) <- unit_ids} else
    if (is.matrix(x)) {rownames(x) <- unit_ids} else
      if (is.data.frame(x)) {row.names(x) <- unit_ids} else {
        stop(
          paste0(
            "Y,W, and X must be vectors, matrices, or perhaps data.frames. ",
            "Factors should be avoided."
          )
        )
      }
  x
}
