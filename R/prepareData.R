
prepareData <- function(
  Y, W, X, match_on, trimming=NULL, model_options
){

  #probably move this into prepareData
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

  if (trimming==0) {analysis_idx <- NULL} else
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
        temp <- utils::capture.output(PF.out <- multinom(W.ref~X))
        ## TODO: make this into a subfunction?
        PF.fit <- stats::fitted(PF.out)

        ## identify sufficient overlap
        overlap_idx <- overlap(PF.fit)$idx
        W <- W[overlap_idx]
        X <- X[overlap_idx,,drop=FALSE]
        Y <- Y[overlap_idx]

        ## organize indices into 'kept' and 'dropped' for clean output
        all_indices <- 1:length(Y)
        analysis_idx <- list(
          indices_kept = overlap_idx,
          indices_dropped = all_indices[!all_indices%in% overlap_idx]
        )
      }
    } else {stop("trimming must be set to FALSE or TRUE")}



  ## Order the treatment increasingly
  ordered_data_list <- reorderByTreatment(W=W,X=X,Y=Y)

  if ( match_on == "existing" ) {
    ## Check that "existing" X has the correctly specified number of levels

    warning("user-supplied propensity scores has not been passing checks as of 2017-03-26")
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
