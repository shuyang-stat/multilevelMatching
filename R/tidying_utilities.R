
#' order the treatment increasingly
#'
#' @param W a treatment vector (1 x n) with numerical values indicating treatment groups
#' @param X a covariate matrix (p x n) with no intercept
#' @param Y a continuous response vector (1 x n)
#'
#' @return the following three elements, ordered according to levels of W
#' \itemize{
#'
#'  \item W: a treatment vector (1 x n) with numerical values indicating treatment groups
#'
#'  \item X: a covariate matrix (p x n) with no intercept
#'
#'  \item Y: a continuous response vector (1 x n)
#
#' }
reorderByTreatment <- function(Y,W,X,method){

  if (method = "MatchOnX") {
    X <- as.matrix(X)
  }
  if (1-is.unsorted(W)) {
    temp <- sort(W,index.return=TRUE)
    temp <- list(x=temp)
    temp$ix <- 1:length(W)
  }
  if (is.unsorted(W)) {
    temp <- sort(W,index.return=TRUE)
  }

  # temp <- orderTrt(W)
  W <- W[temp$ix]
  X <- X[temp$ix,]
  Y <- Y[temp$ix]

  list_to_return <- list(
    W = W,
    X = X,
    Y = Y
  )
  return(list_to_return)
}


argChecks <- function(Y,W,X,method, N=NULL) {



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


#' Strcutres main output in a "tidy" datframe
#'
#' @param untidy_output A list of objects to output from
#' \code{\link{multilevelXMatch}}
#' \code{\link{multilevelGPSMatch}}
#' \code{\link{multilevelGPSStratification}}
#'
#'
#'@return a list containing
#' \itemize{
#'
#'  \item results: A tidy dataframe continaing the following columns:
#'    \itemize{
#'
#'     \item tauestimate:  a vector of estimates for pairwise treatment effects
#'
#'     \item varestimate:  a vector of variance estimates for tauestimate,
#'     using Abadie & Imbens (2006)'s method
#'
#'     \item varestimateAI2012:  a vector of variance estimates for
#'     tauestimate, when matching on the generalized propensity score, using
#'     Abadie & Imbens (2012)'s method.  This variance estimate takes into account
#'     of the uncertainty in estimating the GPS. This will be returned only for
#'     the \code{\link{multilevelGPSMatch}} function and fitting propensity
#'     scores with multinomial logistic regression.
#'
#'    }
#'  \item analysisidx: the index of units after trimming based on Crump et
#'  al. (2009)'s method. This will be returned only if trimming is used in
#'  \code{\link{multilevelGPSMatch}} or \code{\link{multilevelGPSStratification}}
#'
#' }
#'
#'
tidyOutput <- function(untidy_output){

  tidy_dfm <- data.frame(
    stringsAsFactors=FALSE,
    row.names = NULL,
    Param = names(untidy_output$tauestimate),
    Estimate = untidy_output$tauestimate,
    Variance = untidy_output$varestimate
  )

  if ("varestimateAI2012" %in% names(untidy_output)) {
    if (!all(is.na(untidy_output$varestimateAI2012))) {
      tidy_dfm$varestimateAI2012 <- untidy_output$varestimateAI2012
    }
  }

  list_to_return <- list(
    results = tidy_dfm
  )

  if ("analysisidx" %in% names(untidy_output)) {
    list_to_return$analysisidx <- untidy_output$analysisidx
  }

  return(list_to_return)
}
