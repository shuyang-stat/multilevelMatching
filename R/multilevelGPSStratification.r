#' Stratification on GPS with multilevel treatments
#'
#' @inheritParams multilevelMatchX
#' @param GPSM An indicator of the methods used for estimating GPS, options
#'   include "multinomiallogisticReg", "ordinallogisticReg", and "existing"
#' @param NS The number of strata: (only required in the function
#'   \code{\link{multilevelGPSStratification}})
#' @param linearp An indicator of subclassification on GPS (=0) or linear
#'   predictor of GPS (=1): (only required in the function
#'   \code{\link{multilevelGPSStratification}})
#' @param nboot The number of boot replicates for variance estimation: (only
#'   required in the function \code{\link{multilevelGPSStratification}})
#'
#' @return A list with two elements,
#'   \code{tauestimate}, \code{varestimate}, where \code{tauestimate} is a
#'   vector of estimates for pairwise treatment effects, and \code{varestimate}
#'   is a vector of variance estimates, using bootstrapping method.
#'
#' @references Yang, S., Imbens G. W., Cui, Z., Faries, D. E., & Kadziola, Z.
#'   (2016) Propensity Score Matching and Subclassification in Observational
#'   Studies with Multi-Level Treatments. Biometrics, 72, 1055-1065.
#'   \url{https://doi.org/10.1111/biom.12505}
#'
#'   Abadie, A., & Imbens, G. W. (2006). Large sample properties of matching
#'   estimators for average treatment effects. Econometrica, 74(1), 235-267.
#'   \url{https://doi.org/10.1111/j.1468-0262.2006.00655.x}
#'
#'   Abadie, A., & Imbens, G. W. (2016). Matching on the estimated propensity
#'   score. Econometrica, 84(2), 781-807.
#'   \url{https://doi.org/10.3982/ECTA11293}
#'
#' @seealso \code{\link{multilevelGPSMatch}}; \code{\link{multilevelMatchX}}
#'
#' @examples
#'
#' simulated_data <- multilevelMatching::simulated_data
#' set.seed(123)
#' multilevelGPSStratification(
#'   Y = simulated_data$outcome ,
#'   W = simulated_data$treatment,
#'   X = simulated_data[ ,names(simulated_data) %in% paste0("covar", 1:6)],
#'   GPSM = "multinomiallogisticReg",
#'   NS = 5,
#'   linearp = TRUE,
#'   nboot = 10
#' )
#'
#' @export
multilevelGPSStratification <- function(
  Y,W,X,NS,GPSM="multinomiallogisticReg",linearp=0,nboot
  ){

  N <- length(Y) # number of observations
  X <- as.matrix(X)

  trtnumber <- length(unique(W)) # number of treatment levels
  trtlevels <- unique(W) # all treatment levels
  pertrtlevelnumber <- table(W) # number of observations by treatment level
  taunumber <- trtnumber*(trtnumber+1)/2-trtnumber  # number of pairwise treatment effects

  #GPS modeling
  if(GPSM=="multinomiallogisticReg"){
    W.ref <- stats::relevel(as.factor(W),ref="1")
    temp <- utils::capture.output(PF.out <- nnet::multinom(W.ref~X))
    PF.fit <- stats::fitted(PF.out)
    if(linearp==1){
      beta <- stats::coef(PF.out)
      Xbeta <- X%*%t(beta[,-1])
      PF.fit[,-1] <- Xbeta
    }
  }
  if (GPSM == "ordinallogisticReg") {
    PF.out <- MASS::polr(as.factor(W)~X)
    PF.fit <- stats::fitted(PF.out)
  }
  if (GPSM == "existing") {
    PF.fit <- X
  }

  meanwj<-numberwj<-matrix(NA,trtnumber,NS)

  for(kk in 1:trtnumber){
    pwx<-PF.fit[,kk]
    ranking <- stats::ave(pwx,FUN=function(x){
      cut(x,stats::quantile(pwx,(0:NS)/NS,type=2),include.lowest=TRUE,right=FALSE)})
    #type=2 to have the same quintiles as in SAS
    #right=FALSE to have left side closed intervals
    for(jj in 1:NS){
      id<-which(W==kk&ranking==jj)
      meanwj[kk,jj]<-mean(Y[id])
      numberwj[kk,jj]<-sum(ranking==jj) # to get weighted sum at the end
    }
    numberwj<-numberwj*(1-is.na(meanwj))
    meanw<-apply(meanwj*numberwj,1,sum,na.rm=TRUE)/apply(numberwj,1,sum) # to get weighted sum at the end
  }

  tauestimate<-rep(NA,taunumber)

  cnt<-0
  cname1<-c()
  for(jj in 1:(trtnumber-1)){
    for(kk in (jj+1):trtnumber){
      cnt<-cnt+1
      thistrt<-trtlevels[jj]
      thattrt<-trtlevels[kk]
      cname1<-c(cname1,paste(paste(
        paste(paste(paste("EY(",thattrt,sep=""),")",sep=""),
              "-EY(",sep=""),thistrt,sep=""),")",sep=""))
      tauestimate[cnt]<-meanw[kk]-meanw[jj]
    }
  }
  names(tauestimate)<-cname1


  ## Bootstrap the variance
  data<-cbind(W,Y,X)
  results <- boot::boot(data=data, statistic=estforboot,R=nboot,
                  GPSM=GPSM,linearp=linearp,trtnumber=trtnumber,
                  trtlevels=trtlevels,taunumber=taunumber,NS=NS)
  bootvar <- apply(results$t,2,stats::var,na.rm = TRUE)
  names(bootvar)<-cname1

  return(list(tauestimate=tauestimate,varestimate=bootvar))
}
