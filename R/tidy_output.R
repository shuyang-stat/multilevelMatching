


#' DEPRECATED.Strcutres main output in a "tidy" datframe
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
#'  \item analysis_idx: the index of units after trimming based on Crump et
#'  al. (2009)'s method. This will be returned only if trimming is used in
#'  \code{\link{multilevelGPSMatch}} or \code{\link{multilevelGPSStratification}}
#'
#' }
#'
tidyOutput <- function(untidy_output, analysis_idx=NULL){

  stop("tidyOutput is a deprecated utility function and will be removed soon")
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

  # if ("analysis_idx" %in% names(untidy_output)) {
  #   list_to_return$analysis_idx <- untidy_output$analysis_idx
  # }
  #
  if (!is.null(analysis_idx)) {lits_to_return$analysis_idx <- analysis_idx}

  return(list_to_return)
}


nameContrast <- function(trt1,trt2){ paste0("EY(", trt2,")-EY(", trt1,")") }
