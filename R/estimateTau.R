#' Estimate the tau contrasts, perhaps with variance estimates
#'
#' This is a major plumbing function for the package
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
#'
#' @seealso \code{\link{multiMatch}};
#'
#' @return a list including the dataframes for estimates for tau and for mu
#'
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
        Yiw[,kk],
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

# #' Computes Estimated Asymptotic Variance
# #'
# #' See Theorem 7 of Abadie and Imbens 2006 Econometrica
estVarAI2006 <- function(
  N, W, M_matches,
  trt_level_1, trt_level_2,
  Yiw1, Yiw2, tau, Kiw, sigsqiw
){

  Y_contrasts <- Yiw1-Yiw2-tau
  K_M_factor <- (Kiw/M_matches)^2 +
    ((2*M_matches-1)/(M_matches)) *
    ((Kiw/M_matches))
  W_indicator <- (W == trt_level_1 | W == trt_level_2)
  ## From Theorem 7, page 251 of Abadie and Imbens 2006 Econometrica
  V_hat <- mean( Y_contrasts^2 ) + mean( K_M_factor * sigsqiw * W_indicator )


  estimated_asymptotic_variance <- (1/N)*(V_hat)
  estimated_asymptotic_variance
}
