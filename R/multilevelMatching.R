#' Propensity Score Matching and Subclassification in Observational Studies with
#' Multi-level Treatments
#'
#' \pkg{multilevelMatching} implements the estimators introduced in Yang et al.
#' (2016) \emph{Propensity Score Matching and Subclassification in Observational
#' studies with Multi-level Treatments}:
#' \url{https://doi.org/10.1111/biom.12505}. These are covariate- and propensity
#' score-matching estimators for estimating the causal effect of multilevel
#' treatment (i.e., 3 or more treatment types).
#'
#' The main function for estimation via matching on covariates or propensity
#' scores is \code{\link{multiMatch}}. To carry out estimation via
#' subclassification, use \code{\link{multilevelGPSStratification}}.
#'
#'
#' @examples
#'   sim_data <- multilevelMatching::simulated_data
#'   Y <- sim_data$outcome
#'   W <- sim_data$treatment
#'   X <- as.matrix(sim_data[ ,-(1:2)])
#'   names(Y) <- paste0("ID", 1:length(Y))
#'
#'   trimming <- FALSE
#'   method <- c("covariates", "polr", "multinom")[2]
#'
#'   multiMatch(Y,W,X,trimming=trimming,match_on=method)
#'
#'
#' @name multilevelMatching
#' @docType package
NULL



