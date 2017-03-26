#' Stratification on GPS with multilevel treatments
#'
#' @param Y a continuous response vector (1 x n)
#' @param W a treatment vector (1 x n) with numerical values indicating treatment groups
#' @param X a covariate matrix (p x n) with no intercept
#' @param GPSM an indicator of the methods used for estimating GPS, options include "multinomiallogisticReg", "ordinallogisticReg", and "existing"
#' @param NS (only required in the function multilevelGPSStratification) the number of strata
#' @param linearp (only required in the function multilevelGPSStratification) an indicator of subclassification on GPS (=0) or linear predictor of GPS (=1)
#' @param nboot (only required in the function multilevelGPSStratification) the number of boot replicates for variance estimation
#'
#' @return according to \code{\link{tidyOutput}}: a dataframe with two columns,
#' tauestimate, varestimate, where
#' tauestimate is a vector of estimates for pairwise treatment effects, and
#' varestimate is a vector of variance estimates, using bootstrapping method.
#'
#' @seealso \code{\link{multilevelGPSMatch}}; \code{\link{multilevelMatchX}}
#'
#'
#' @import Matching boot nnet optmatch MASS
#'
#' @export
multilevelGPSStratification <- function(Y,W,X,NS,GPSM="multinomiallogisticReg",linearp=0,nboot){
  ### Generalized propensity score stratification
  ## PF.out=Propensity function estimation
  ## linearp=1 if use linear predictor for stratification
  ## X=design matrix
  ## Y=outcome
  ## W=treatment indicator
  ## nboot: the number of bootstrap replicates in variance estimator
  ## return values: point estimator,bootstrapping variance estimator
  #

  ## some checks
  match_method <- "StratifyOnGPS"
  argChecks(Y=Y, W=W, X=X, match_method = match_method, N=NULL)

  N <- length(Y) # number of observations
  ## order the treatment increasingly
  ordered_data <- reorderByTreatment(W=W,X=X,Y=Y)
  W <- ordered_data$W
  X <- ordered_data$X
  Y <- ordered_data$Y
  ## some checks, again
  argChecks(Y=Y, W=W, X=X, match_method = match_method, N=N)


  trtnumber<-length(unique(W)) # number of treatment levels
  trtlevels<-unique(W) # all treatment levels
  pertrtlevelnumber<-table(W) # number of observations by treatment level
  taunumber<-trtnumber*(trtnumber+1)/2-trtnumber  # number of pairwise treatment effects

  #PF modeling
  if(GPSM=="multinomiallogisticReg"){
    W.ref <- relevel(as.factor(W),ref="1")
    temp<-capture.output(PF.out <- multinom(W.ref~X))
    PF.fit <- fitted(PF.out)
    if(linearp==1){
      beta<-coef(PF.out)
      Xbeta<-X%*%t(beta[,-1])
      PF.fit[,-1]<-Xbeta
    }
  }
  if(GPSM=="ordinallogisticReg"){
    PF.out <- polr(as.factor(W)~X)
    PF.fit <- fitted(PF.out)
  }
  if(GPSM=="existing"){
    if (any(dim(X)!=c(1,N))){
      stop("user-supplied propensity scores (through argument 'X') should
           be a vector of same length as Y and W")
    }
    #need to check the row sum of X is 1 - debug1
    PF.fit <- X
  }

  meanwj<-numberwj<-matrix(NA,trtnumber,NS)

  for(kk in 1:trtnumber){
    pwx<-PF.fit[,kk]
    ranking<-ave(pwx,FUN=function(x)cut(x,quantile(pwx,(0:NS)/NS,type=2),include.lowest=TRUE,right=FALSE))
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
      cname1<-c(cname1,paste(paste(paste(paste(paste("EY(",thattrt,sep=""),")",sep=""),"-EY(",sep=""),thistrt,sep=""),")",sep=""))
      tauestimate[cnt]<-meanw[kk]-meanw[jj]
    }
  }
  names(tauestimate)<-cname1


  data<-cbind(W,Y,X)
  results <- boot(data=data, statistic=estforboot,R=nboot,
                  GPSM=GPSM,linearp=linearp,trtnumber=trtnumber,
                  trtlevels=trtlevels,taunumber=taunumber,NS=NS)
  bootvar<-apply(results$t,2,var,na.rm = TRUE)
  names(bootvar)<-cname1

  untidy_output <- list(tauestimate=tauestimate,varestimate=bootvar)
  tidy_output <- tidyOutput(untidy_output=untidy_output)
  return(tidy_output)
}
