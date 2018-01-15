
# #' for naming columns in the matching matrix
nameCols <- function(trtlevels){
  col_names <- lapply(trtlevels, function(x){
    c(paste0("m",x, ".1"), paste0("m", x, ".2"))
  })
  unlist(col_names)
}



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




# #' order the treatment increasingly
# #'
# #' @param W a treatment vector (1 x n) with numerical values indicating treatment groups
# #' @param X a covariate matrix (p x n) with no intercept
# #' @param Y a continuous response vector (1 x n)
# #'
# #' @return the following elements, ordered according to levels of W
# #' \itemize{
# #'
# #'  \item W: a treatment vector (1 x n) with numerical values indicating treatment groups
# #'
# #'  \item X: a covariate matrix (p x n) with no intercept
# #'
# #'  \item Y: a continuous response vector (1 x n)
# #
# #' }
# #' along with these downstream elements of treatment:
# #' \itemize{
# #' \item trtnumber: number of treatment levels
# #' \item trtlevels: all treatment levels
# #' \item pertrtlevelnumber: number of observations by treatment level
# #' \item taunumber: number of pairwise treatment effects
# #' }
reorderByTreatment <- function(Y,W,X){

  N <- length(Y)

  if (1-is.unsorted(W)) {
    temp <- sort(W, index.return=TRUE)
    temp <- list(x = temp)
    ## TODO: recode this into
    ## temp <- list(ix = 1:length(W))
    ## (add unit test before attempting)
    temp$ix <- 1:length(W)
  }
  if (is.unsorted(W)) {
    temp <- sort(W,index.return=TRUE)
  }

  orig_to_sorted <- temp$ix
  sorted_to_orig <- order(temp$ix)

  W <- W[orig_to_sorted]
  X <- X[orig_to_sorted, , drop=FALSE]
  Y <- Y[orig_to_sorted]

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

  list_to_return
}
