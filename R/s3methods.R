

### Methods for the multiMatch class


#' Prints a summary of the estimates from a multiMatch object
#'
#' @param x object of class "multiMatch"
#' @param ... dots
#'
#' @method print multiMatch
#'
#' @examples
#'
#'   sim_data <- multilevelMatching::simulated_data
#'   Y <- sim_data$outcome
#'   W <- sim_data$treatment
#'   X <- as.matrix(sim_data[ ,-(1:2)])
#'   names(Y) <- paste0("ID", 1:length(Y))
#'
#'   trimming <- FALSE
#'   method <- c("covariates", "polr", "multinom")[2]
#'
#'   fit <- multiMatch(Y,W,X,trimming=trimming,match_on=method)
#'   print(fit)
#'
#' @export
print.multiMatch <- function(x,...){

  ests <- x$results
  # idx <- x$analysis_idx
  estimate_args <- x$estimate_args

  cat("-------------- Causal estimates ---------------\n")

  if (all(is.na(ests$VarianceAI2016))){
    ests$VarianceAI2016 <- NULL
  }
  print(ests)

  second_message <- paste0(
    "--- Matching on '", estimate_args$match_on,
    "' with M=",
    estimate_args$M_matches,
    ", J=",
    estimate_args$J_var_matches,
    " ---\n"
  )
  cat(second_message)
}


#' Prints a summary of a multiMatch object
#'
#' @param object object of class "multiMatch"
#' @param ... dots
#'
#' @method summary multiMatch
#'
#' @author Brian G. Barkley
#'
#'
#' @examples
#'
#'   sim_data <- multilevelMatching::simulated_data
#'   Y <- sim_data$outcome
#'   W <- sim_data$treatment
#'   X <- as.matrix(sim_data[ ,-(1:2)])
#'   names(Y) <- paste0("ID", 1:length(Y))
#'
#'   trimming <- FALSE
#'   method <- c("covariates", "polr", "multinom")[2]
#'
#'   fit <- multiMatch(Y,W,X,trimming=trimming,match_on=method)
#'   summary(fit)
#'
#' @export
summary.multiMatch <- function(object, ...){


  ests <- object$results
  # idx <- object$analysis_idx
  estimate_args <- object$estimate_args

  cat("------------- Method arguments --------------\n")

  args2print <- estimate_args[c(
    "match_on", "model_options", "M_matches", "J_var_matches",
    "trt_levels", "N_per_trt"
  )]

  print(args2print)

  cat("------------- Causal estimates --------------\n")


  if (all(is.na(ests$VarianceAI2016))){
    ests$VarianceAI2016 <- NULL
  }
  print(ests)

  cat("---------------------------------------------\n")
}




### More methods for the legacy functions


