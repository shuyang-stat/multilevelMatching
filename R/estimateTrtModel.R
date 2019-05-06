
#' Estimate Treatment Model for Generalized Propensity Scores
#'
#' This function is used to fit the model for the generalized propensity score.
#' Users can apply this function before \code{\link{multiMatch}} and verify that
#' the output's fitted model object is the same as the user desires.
#'
#' Note that the \code{model_options} argument must be a list with
#' \code{reference_level} element. Future versions of this package may allow
#' for the user to supply a fitted model object directly to
#' \code{\link{multiMatch}}; to request this feature, users should go to the
#' GitHub repository and fill out an Issue requesting it.
#'
#' @inheritParams multiMatch
#' @param ... the dots argument
#'
#' @return A list element with two items: \itemize{
#'   \item \code{prop_score_model} the fitted model object
#'   \item \code{prop_score_ests} the estimated generalized propensity scores
#'   for each individual in the dataset
#' }
#'
#' @export
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
#'
#' trt_model <- do.call(estimateTrtModel, prepared_data)
#' estimated_generalized_propensity_scores <- trt_model$prop_score_ests
#'
estimateTrtModel <- function(
  W, X, match_on, model_options, ...
){
  if (match_on == "multinom") {
    W.ref <- stats::relevel(as.factor(W),ref=model_options$reference_level)
    temp <- utils::capture.output(prop_score_model <- nnet::multinom(W.ref~X))
    prop_score_ests <- stats::fitted(prop_score_model)
    vcov_coeff <- stats::vcov(prop_score_model)
  } else
    if (match_on == "polr") {
      prop_score_model <- MASS::polr(as.factor(W)~X)
      prop_score_ests <- stats::fitted(prop_score_model)
    } else
      if (match_on == "existing") {
        prop_score_model <- NULL
        prop_score_ests <- X
      } else
        if (match_on == "covariates"){
          prop_score_model <- NULL
          prop_score_ests <- NULL
        } else {
          stop("match_on not recognized")
        }
  out_list <- list(
    prop_score_model = prop_score_model,
    prop_score_ests  = prop_score_ests
  )
  if (match_on == "multinom"){out_list$vcov_coeff <- vcov_coeff}
  out_list
}
