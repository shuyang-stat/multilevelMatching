#' Simulated dataset for multilevelMatching package
#'
#' This is a dataset with six baseline covariates, one column indicating
#' treatment level, and one column indicating post-treatment outcome. This
#' simulated data is purely for illustration, and any combination of the
#' covariates can be assumed to sufficient to meet conditional exchangeability.
#'
#' @format A data frame with 300 rows and 8 variables:
#'   \describe{
#'     \item{outcome}{Outcome of interest}
#'     \item{treatment}{Treatment level of the unit}
#'     \item{covar1}{Baseline covariate 1}
#'     \item{covar2}{Baseline covariate 2}
#'     \item{covar3}{Baseline covariate 3}
#'     \item{covar4}{Baseline covariate 4}
#'     \item{covar5}{Baseline covariate 5}
#'     \item{covar6}{Baseline covariate 6}
#'  }
#' @source \url{http://www.diamondse.info/}
"simulated_data"
