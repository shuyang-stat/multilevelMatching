#' Estimate the tau contrasts, perhaps with variance estimates
#'
#' This is a major plumbing function for the package
#'
#' @param trtlevels vector of the unique levels of treatment W
#' @param meanw vector of the estimated mean w.r.t. each treatment w
#' @param trtnumber a scalar, the number of treatment levels
#' @param taunumber a scalar, the number of tau contrasts to estimate
#' @param N A scalar for the number of rows in the data
#' @param Yiw for AI2012 variance
#' @param Kiw for AI2012 variance
#' @param sigsqiw for AI2012 variance
#' @param W for AI2012 variance
#'
#' To estimate variance in Matching using AI2006 method in MatchX and MatchGPS,
#' all of the following must be supplied. Otherwise, they default to NULL.
#'   \itemize{
#'
#'   \item @param Yiw:  a vector of estimates for pairwise treatment effects
#'
#'   \item @param Kiw
#'
#'   \item @param sigsqiw
#'
#'   \item @param W
#'   }
#'
#'  @seealso \code{\link{multilevelGPSMatch}}; \code{\link{multilevelMatchX}};
#'  \code{\link{multilevelGPSStratification}}
#'
#' @return a list including the dataframes for estimates for tau and for mu
#'
#'
#' @export
estimateTau <- function(
  trtlevels,meanw,
  trtnumber,taunumber,N,
  Yiw=NULL, Kiw=NULL,sigsqiw=NULL,W=NULL
){

  tau_dfm <- list(
    # stringsAsFactors = FALSE, row.names = NULL,
    Param = rep(NA, taunumber),
    Trt1 = rep(NA, taunumber),
    Trt2 = rep(NA, taunumber),
    Estimate = rep(NA, taunumber),
    Variance = rep(NA, taunumber),
    VarianceAI2012 = rep(NA,taunumber)
  )


  mu_dfm <- data.frame(
    stringsAsFactors = FALSE, row.names = NULL,
    Param = nameMu(trtlevels),
    Trt = trtlevels,
    Estimate = meanw
  )
  imputes_mat <- Yiw#data.frame(
  #   stringsAsFactors = FALSE, row.names = NULL,
  #   Param = nameMu(trtlevels),
  #   Trt = trtlevels,
  #   Estimate = meanw
  # )

  row_num <- 0
  # cname1<-c()
  for(jj in 1:(trtnumber-1)){
    for(kk in (jj+1):trtnumber){
      row_num <- row_num+1
      # thistrt <- trtlevels[jj]
      # thattrt <- trtlevels[kk]

      tau_dfm$Trt1[row_num] <- trtlevels[jj]
      tau_dfm$Trt2[row_num] <- trtlevels[kk]
      tau_dfm$Param[row_num] <- nameContrast(trt1=tau_dfm$Trt1[row_num], trt2=tau_dfm$Trt2[row_num])
      tau_dfm$Estimate[row_num] <- meanw[kk]-meanw[jj]
      tau_dfm$Variance[row_num] <- (1/N)*(
        mean( (Yiw[,kk]-Yiw[,jj]-(tau_dfm$Estimate[row_num]))^2 ) +
          mean( (Kiw^2+Kiw)*sigsqiw*
                  (W==tau_dfm$Trt1[row_num] | W==tau_dfm$Trt2[row_num])
          )
      )

      # varestimate[row_num]<-mean((Yiw[,kk]-Yiw[,jj]-(meanw[kk]-meanw[jj]))^2)+
      # mean((Kiw^2+Kiw)*sigsqiw*(W==thistrt | W==thattrt))
    }
  }
  if (row_num != taunumber) { stop("Error in for loop") }

  tau_dfm <-
    as.data.frame(tau_dfm, stringsAsFactors = FALSE, row.names = NULL)
  results <- list(
    tau_dfm = tau_dfm,
    mu_dfm = mu_dfm
  )
  return(results)
}
