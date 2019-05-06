
# #' Name Two or All Columns in the Matching Matrix
# #'
# #' Useful when using Abadie and Imbens (2016) variance estimator
# #'
# #' @inheritParams matchAllTreatments
nameCols <- function(trt_levels){
  col_names <- lapply(trt_levels, function(x){
    c(paste0("m",x, ".1"), paste0("m", x, ".2"))
  })
  unlist(col_names)
}



#' Naming the matching contrasts
#'
#' @param trt1 Former treatment level
#' @param trt2 Latter treatment level
#'
#' @examples
#'   nameContrast(trt1=1, trt2=0)
#'
#' @export
nameContrast <- function(trt1,trt2){ paste0("EY(", trt2,")-EY(", trt1,")") }

#' Naming the matching population mean mu's
#'
#' @param trt Treatment level
#'
#' @examples
#'  nameMu(1)
#'
#' @export
nameMu <- function(trt){ paste0("EY(", trt,")") }


# #' Defensive programming for data re-ordering
# #'
# #' @inheritParams multiMatch
# #' @param N The number of unique units
argChecks <- function(Y,W,X,N=NULL) {

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




# #' Orders the treatment increasingly
# #'
# #' @inheritParams setIDs
# #' @param W A treatment vector (1 x n) with numerical values indicating treatment groups
# #' @param X A covariate matrix (p x n) with no intercept
# #' @param Y A continuous response vector (1 x n)
# #' @param unit_ids_unsorted The \code{unit_ids} before the data is reordered
# #'
# #' @return The following elements, ordered according to levels of \code{W}
# #' \itemize{
# #'
# #'  \item \code{W}: a treatment vector (1 x n) with numerical values indicating treatment groups
# #'
# #'  \item \code{X}: a covariate matrix (p x n) with no intercept
# #'
# #'  \item \code{Y}: a continuous response vector (1 x n)
# #
# #' }
# #' along with these downstream elements of treatment:
# #' \itemize{
# #' \item \code{num_trts}: number of treatment levels
# #' \item \code{trt_levels}: all treatment levels
# #' \item \code{N_per_trt}: number of observations by treatment level
# #' \item \code{num_contrasts}: number of pairwise treatment effects
# #' \item \code{orig_to_sorted}: vector to rearrange from original to sorted by treatment
# #' \item \code{sorted_to_orig}: vector to rearrange from sorted to original order
# #' }
reorderByTreatment <- function(Y,W,X, unit_ids_unsorted){

  N <- length(Y)

  if (!is.unsorted(W)) {
    # temp <- sort(W, index.return=TRUE)
    ## Above does not work with factor variables
    ## (See test_utilities.R)
    temp <- list(x = W)
    temp$ix <- 1:length(W)
  } else
  if (is.unsorted(W)) {
    temp <- sort(W,index.return=TRUE)
  }

  orig_to_sorted <- temp$ix
  sorted_to_orig <- order(temp$ix)

  ## Defensive programming
  stopifnot(identical(getIDs(W), unit_ids_unsorted))
  stopifnot(identical(getIDs(Y), unit_ids_unsorted))
  stopifnot(identical(getIDs(X), unit_ids_unsorted))

  W <- W[orig_to_sorted]
  Y <- Y[orig_to_sorted]
  X <- X[orig_to_sorted, , drop=FALSE]


  unit_ids_sorted <- unit_ids_unsorted[orig_to_sorted]

  W <- setIDs(W, unit_ids_sorted)
  Y <- setIDs(Y, unit_ids_sorted)
  X <- setIDs(X, unit_ids_sorted)

  ## Defensive programming
  stopifnot(identical(getIDs(W)[sorted_to_orig], unit_ids_unsorted))
  stopifnot(identical(getIDs(Y)[sorted_to_orig], unit_ids_unsorted))
  stopifnot(identical(getIDs(X)[sorted_to_orig], unit_ids_unsorted))

  ## Defensive checks, again
  argChecks(Y = Y, W = W, X = X, N = N)

  list_to_return <- list(W = W, X = X, Y = Y, N = N)

  ## Adding the following to output
  num_trts <- length(unique(W))
  list_to_return$num_trts <- num_trts # number of treatment levels
  list_to_return$trt_levels <- unique(W) # all treatment levels
  list_to_return$N_per_trt <- table(W) # number of observations by treatment level
  list_to_return$num_contrasts <- (num_trts*(num_trts+1)/2)-num_trts  # number of pairwise treatment effects
  list_to_return$orig_to_sorted <- orig_to_sorted
  list_to_return$sorted_to_orig <- sorted_to_orig
  list_to_return$unit_ids_unsorted <- unit_ids_unsorted
  list_to_return$unit_ids_sorted <- unit_ids_sorted


  list_to_return
}
