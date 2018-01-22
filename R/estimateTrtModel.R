
#' Estimate Treatment Model for Generalized Propensity Scores
#'
#' Note that the model_options argument must be a list with reference_level
#' element
#'
#' @inheritParams multiMatch
#' @param ... the dots argument
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
