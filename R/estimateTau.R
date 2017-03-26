#' Estimate the tau contrasts, perhaps with variance estimates
#'
#' This is a major plumbing function for the package
#'
#' @param trtlevels vector of the unique levels of treatment W
#' @param meanw vector of the estimated mean w.r.t. each treatment w
#' @param trtnumber a scalar, the number of treatment levels
#' @param taunumber a scalar, the number of tau contrasts to estimate
#' @param N A scalar for the number of rows in the data
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
#' @return a tidy dataframe including Estimate, Variance, and referencing the
#'   target parameter (tau)
#'
#' @export
estimateTau <- function(
  trtlevels,meanw,
  trtnumber,taunumber,N,
  Yiw=NULL, Kiw=NULL,sigsqiw=NULL,W=NULL
){

  results_dfm <- list(
    Param = rep(NA, taunumber),
    Trt1 = rep(NA, taunumber),
    Trt2 = rep(NA, taunumber),
    Estimate = rep(NA, taunumber),
    Variance = rep(NA, taunumber),
    VarianceAI2012 = rep(NA,taunumber)
  )

  row_num <- 0
  # cname1<-c()
  for(jj in 1:(trtnumber-1)){
    for(kk in (jj+1):trtnumber){
      row_num <- row_num+1
      # thistrt <- trtlevels[jj]
      # thattrt <- trtlevels[kk]

      results_dfm$Trt1[row_num] <- trtlevels[jj]
      results_dfm$Trt2[row_num] <- trtlevels[kk]
      results_dfm$Param[row_num] <- nameContrast(trt1=results_dfm$Trt1[row_num], trt2=results_dfm$Trt2[row_num])
      results_dfm$Estimate[row_num] <- meanw[kk]-meanw[jj]
      results_dfm$Variance[row_num] <- (1/N)*(
        mean( (Yiw[,kk]-Yiw[,jj]-(results_dfm$Estimate[row_num]))^2 ) +
          mean( (Kiw^2+Kiw)*sigsqiw*
                  (W==results_dfm$Trt1[row_num] | W==results_dfm$Trt2[row_num])
          )
      )

      # varestimate[row_num]<-mean((Yiw[,kk]-Yiw[,jj]-(meanw[kk]-meanw[jj]))^2)+
      # mean((Kiw^2+Kiw)*sigsqiw*(W==thistrt | W==thattrt))
    }
  }
  if (row_num != taunumber) { stop("Error in for loop") }

  results_dfm <-
    as.data.frame(results_dfm, stringsAsFactors = FALSE, row.names = NULL)
  return(results_dfm)
}
