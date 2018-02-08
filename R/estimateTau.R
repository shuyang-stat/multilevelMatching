#' Calculate the estimates of population-level estimands (e.g., tau).
#'
#' This is a major plumbing function for the package. All matching procedures
#' are carried out in \code{\link{matchImputePO}} (for point estimates) and
#' \code{\link{estSigSq}} (for variance), which are subfunctions of
#' \code{\link{matchAllTreatments}}. Most of the necessary arguments to this
#' function are output from these two subfunctions.
#'
#' @inheritParams estimateTrtModel
#' @inheritParams multiMatch
#' @param trt_levels vector of the unique levels of treatment W
#' @param num_trts a scalar, the number of treatment levels
#' @param num_contrasts a scalar, the number of tau contrasts to estimate
#' @param N A scalar for the number of rows in the data
#' @param Yiw Matrix of all imputed potential outcomes
#' @param mean_Yiw vector of the estimated mean w.r.t. each treatment w
#' @param Kiw Vector of times each unit is matched to
#' @param sigsqiw Estimated sigma squared, from AI2006
#'
#' @seealso \code{\link{multiMatch}}
#'
#' @return A list, including the tidy dataframes estimates of target estimands
#'
estimateTau <- function(
  trt_levels, mean_Yiw,
  num_trts, num_contrasts, N, M_matches,
  Yiw, Kiw, sigsqiw, W,
  ...
){

  blank_vec <- rep(NA, num_contrasts)
  tau_dfm <- list(
    # stringsAsFactors = FALSE, row.names = NULL,
    # Easier to construct this as a list object
    Param = blank_vec,
    Trt1 = blank_vec,
    Trt2 = blank_vec,
    Estimate = blank_vec,
    Variance = blank_vec,
    VarianceAI2012 = blank_vec
  )


  mu_dfm <- data.frame(
    Param = nameMu(trt_levels),
    Trt = trt_levels,
    Estimate = mean_Yiw,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  row_num <- 0

  for(jj in 1:(num_trts-1)){
    for(kk in (jj+1):num_trts){
      row_num <- row_num+1

      tau_dfm$Trt1[row_num] <- trt_levels[jj]
      tau_dfm$Trt2[row_num] <- trt_levels[kk]
      tau_dfm$Param[row_num] <- nameContrast(trt1=tau_dfm$Trt1[row_num], trt2=tau_dfm$Trt2[row_num])
      tau_dfm$Estimate[row_num] <- mean_Yiw[kk]-mean_Yiw[jj]
      tau_dfm$Variance[row_num] <- estVarAI2006(
        N = N, W = W, Kiw = Kiw, sigsqiw = sigsqiw, M_matches = M_matches,
        Yiw1 = Yiw[,kk],
        Yiw2 = Yiw[,jj],
        trt_level_1 = tau_dfm$Trt1[row_num],
        trt_level_2 = tau_dfm$Trt2[row_num],
        tau = tau_dfm$Estimate[row_num]
      )
    }
  }

  tau_dfm <-
    as.data.frame(tau_dfm, stringsAsFactors = FALSE, row.names = NULL)
  results <- list(
    tau_dfm = tau_dfm,
    mu_dfm = mu_dfm
  )
  results
}

#' Computes Estimated Asymptotic Variance of matching estimators.
#'
#' See Theorem 7 of Abadie and Imbens 2006 Econometrica for the formula.
#'
#' @inheritParams estimateTau
#' @inheritParams multiMatch
#' @param tau Estimated value of \eqn{\tau(W_1, W_2)}.
#' @param trt_level_1 Unique treatment level 1; aka \eqn{W_1} in \eqn{\tau(W_1,
#'   W_2)}
#' @param trt_level_2 Unique treatment level 2; aka \eqn{W_2} in \eqn{\tau(W_1,
#'   W_2)}
#' @param Yiw1 Vector of imputed outcomes for all units for \code{trt_level_2}.
#' @param Yiw2 Vector of imputed outcomes for all units for \code{trt_level_2}.
#'
#' @return A single numeric value for the estimated asymptotic variance of the
#'   estimator.
#'
estVarAI2006 <- function(
  N, W, M_matches,
  trt_level_1, trt_level_2,
  Yiw1, Yiw2, tau, Kiw, sigsqiw
){

  Y_contrasts <- (Yiw1-Yiw2)-tau
  Y_contrasts_sq <- Y_contrasts^2
  ## Estimating variance of conditional mean
  V_taux <- mean(Y_contrasts_sq)

  K_M_factor <- calcKMFactor(Kiw, M_matches)
  W_indicator <- (W == trt_level_1 | W == trt_level_2)
  ## Estimating conditional variance
  V_E <- mean( K_M_factor * sigsqiw * W_indicator )

  ## Estimating marginal variance
  ## From Theorem 7, page 251 of Abadie and Imbens 2006 Econometrica
  V_hat <- V_taux + V_E
  estimated_asymptotic_variance <- (1/N)*(V_hat)

  estimated_asymptotic_variance
}

#' Calculate the variance component for number of times unit is a match.
#'
#' This function calculates \code{K_M_factor}, a numeric vector. Each entry in
#' this vector is a function of the number of times each unit is matched to, aka
#' \eqn{K_M(i)} (corresponding to \code{Kiw}, where \eqn{M} corresponds to \code{M_matches}. The calculation
#' in this function comes from Theorem 7, page 251 of Abadie and Imbens 2006
#' Econometrica. The \code{K_M_factor} is an important component in the variance
#' estimation, and is called from \code{\link{estVarAI2006}} in
#' \code{\link{estimateTau}}.
#'
#' @inheritParams estVarAI2006
#'
#' @return A numeric vector.
#'
#' This function is exported for use in other packages.
#'
#' @export
calcKMFactor <- function(Kiw, M_matches){
  (Kiw/M_matches)^2 + ( (2*M_matches-1)/(M_matches) ) * (Kiw/M_matches)
}
