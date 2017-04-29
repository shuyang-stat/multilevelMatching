

#' naming the matching contrasts
#' @param trt1 former treatment level
#' @param trt2 latter treatment level
#'
#' @export
nameContrast <- function(trt1,trt2){ paste0("EY(", trt2,")-EY(", trt1,")") }

#' naming the matching mu's
#' @param trt treatment level
#'
#' @export
nameMu <- function(trt){ paste0("EY(", trt,")") }

prepareData <- function(
  Y,W,X,match_method,
  GPSM = NULL,
  Trimming = FALSE#Trimming_fit_args
){

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
  ## some checks
  argChecks(Y=Y, W=W, X=X, N=NULL)


  if (Trimming==0) {analysis_idx <- NULL}#1:length(Y) }

  if (match_method == "MatchOnX") {
    #nothing
  } else if (match_method == "StratifyOnGPS") {
    #nothing
  } else if (match_method == "MatchOnGPS"){

    if(Trimming==1){
      #PF modeling
      W.ref <- relevel(as.factor(W),ref=1)
      temp<-capture.output(PF.out <- multinom(W.ref~X))  ##make this into a subfunction?
      PF.fit <- fitted(PF.out)

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

  } else {stop("match_method coding error")}



  # N <- length(Y) # number of observations ##moved into reorderByTreatment()

  ## order the treatment increasingly
  ordered_data_list <- reorderByTreatment(W=W,X=X,Y=Y)
  ## check that "existing" X has the correctly specified number of levels
  if ((!is.null(GPSM))&&(GPSM=="existing")) {

    warning("user-supplied propensity scores has not been passing checks as of 2017-03-26")
      if (nrow(X) != ordered_data_list$N) {
        stop("user-supplied propensity scores (through argument X and
              GPSM='existing') should be a matrix with number of rows
              equal to length of Y and W")
      }
      if (ncol(X) != ordered_data_list$trtnumber) {
        stop("user-supplied propensity scores (through argument X and
              GPSM='existing') should be a matrix with number of columns
             equal to the number of treatment levels in W")
      }
      #need to check the row sum of X is 1 - debug1
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
      # PF.fit <- X ##this happens later
  }

  ordered_data_list$analysis_idx <- analysis_idx
  return(ordered_data_list)
}



#' order the treatment increasingly
#'
#' @param W a treatment vector (1 x n) with numerical values indicating treatment groups
#' @param X a covariate matrix (p x n) with no intercept
#' @param Y a continuous response vector (1 x n)
#'
#' @return the following elements, ordered according to levels of W
#' \itemize{
#'
#'  \item W: a treatment vector (1 x n) with numerical values indicating treatment groups
#'
#'  \item X: a covariate matrix (p x n) with no intercept
#'
#'  \item Y: a continuous response vector (1 x n)
#
#' }
#' along with these downstream elements of treatment:
#' \itemize{
#' \item trtnumber: number of treatment levels
#' \item trtlevels: all treatment levels
#' \item pertrtlevelnumber: number of observations by treatment level
#' \item taunumber: number of pairwise treatment effects
#' }
reorderByTreatment <- function(Y,W,X){

  N <- length(Y)

  if (1-is.unsorted(W)) {
    temp <- sort(W,index.return=TRUE)
    temp <- list(x=temp)
    ## future todo: recode this into
    ## temp <- list(ix = 1:length(W))
    ## (needs unit test though)
    temp$ix <- 1:length(W)
  }
  if (is.unsorted(W)) {
    temp <- sort(W,index.return=TRUE)
  }

  orig_to_sorted <- temp$ix
  sorted_to_orig <- order(temp$ix)

  # temp <- orderTrt(W)
  W <- W[orig_to_sorted]
  X <- X[orig_to_sorted,,drop=FALSE]
  Y <- Y[orig_to_sorted]

  ## some checks, again
  argChecks(Y=Y,W=W,X=X,N=N)

  list_to_return <- list(W = W,X = X,Y = Y,N=N)

  ## adding to output
  trtnumber <- length(unique(W))
  list_to_return$trtnumber <- trtnumber # number of treatment levels
  list_to_return$trtlevels <- unique(W) # all treatment levels
  list_to_return$pertrtlevelnumber <- table(W) # number of observations by treatment level
  list_to_return$taunumber <- trtnumber*(trtnumber+1)/2-trtnumber  # number of pairwise treatment effects
  list_to_return$orig_to_sorted <- orig_to_sorted
  list_to_return$sorted_to_orig <- sorted_to_orig


  return(list_to_return)
}


argChecks <- function(Y,W,X,#match_method,
                      N=NULL) {



  if ((length(W) != length(Y))) {
    # write a unit test here
    stop("length of Y must equal length of W")
  }
  if ((nrow(X) != length(Y))) {
    # write a unit test here
    stop("length of Y must equal the number of rows in matrix X (or length of X)")
  }
  if ( (!is.null(N)) && length(Y)!=N) {
    stop("Re-ordering data has failed; Y length has changed")
  }
}
