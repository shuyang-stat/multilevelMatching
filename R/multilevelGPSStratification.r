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
#' @return A list with 2 elements: tauestimate, varestimate, where
#' tauestimate is a vector of estimates for pairwise treatment effects, and
#' varestimate is a vector of variance estimates, using bootstrapping method.
#'
#' @seealso \code{\link{multilevelGPSMatch}}; \code{\link{multilevelMatchX}}
#'
#' @examples
#'
#'   set.seed(111)
#'   n    <- 5000*6
#'   # X1-X3 3 MVN var 2, 1, 1, covars 1, -1, -.5
#'   vars   <- c(2,1,1)
#'   covars <- c(1,-1,-.5)
#'   mu     <- c(0,0,0)
#'   tau    <- 1
#'   Sigma <- diag(vars)
#'   Sigma[2,1] <- Sigma[1,2] <- covars[1]
#'   Sigma[3,1] <- Sigma[1,3] <- covars[2]
#'   Sigma[3,2] <- Sigma[2,3] <- covars[3]
#'   trt1 <- 100; trt1
#'   trt2 <- 100; trt2
#'   trt3 <- 100; trt3
#'   # draw Xs
#'   X13 <- mvrnorm(n,mu=mu,Sigma=Sigma, empirical = FALSE)
#'   X1 <- X13[,1]
#'   X2 <- X13[,2]
#'   X3 <- X13[,3]
#'   X4 <- runif(n,-3,3)
#'   X5 <- rchisq(n, df=1)
#'   X6 <- rbinom(n,size=1,prob=.5)
#'
#'   xb2 <- 0.1*(X1^2+X2+X3+X4+X5+X6)
#'   xb3 <- 0.1*(X1+X2^2+X3^2+X4+X5+X6)
#'   exb2<-exp(xb2)
#'   exb3<-exp(xb3)
#'   pi1<-1/(1+exp(xb2)+exp(xb3))
#'   pi2<-exp(xb2)/(1+exp(xb2)+exp(xb3))
#'   pi3<-exp(xb3)/(1+exp(xb2)+exp(xb3))
#'   pi<-cbind(pi1,pi2,pi3)
#'   apply(pi,2,mean)
#'
#'   W<-matrix(NA,n,4)
#'   colnames(W)   <- c("W1","W2","W3","W")
#'   for(kk in 1:n){
#'     W[kk,1:3]<-rmultinom(1, 1, prob = pi[kk,])
#'   }
#'
#'   sim.dat <- data.frame(W,X1,X2,X3,X4,X5,X6)
#'   trt1.keep <- sample(which(sim.dat$W1==1),trt1,replace=FALSE)
#'   trt2.keep <- sample(which(sim.dat$W2==1),trt2,replace=FALSE)
#'   trt3.keep <- sample(which(sim.dat$W3==1),trt3,replace=FALSE)
#'   sim.dat <- sim.dat[c(trt1.keep,trt2.keep,trt3.keep),]
#'   sim.dat[,"W"]<-sim.dat[,"W1"]+2*sim.dat[,"W2"]+3*sim.dat[,"W3"]
#'   sim.dat[,"W"]<-as.factor(sim.dat[,"W"])
#'   W <- sim.dat[,"W"]
#'   X <- as.matrix(sim.dat[,names(sim.dat)[-c(1:4)]])
#'   X1 <- X[,"X1"]; X2 <- X[,"X2"]; X3 <- X[,"X3"]; X4 <- X[,"X4"]; X5 <- X[,"X5"];X6 <- X[,"X6"]
#'
#'   # outcome: treatment effect is zero
#'   u  <- rnorm(nrow(X))
#'   # ouctome 1 (linear)
#'   Y <- 	(W==1)*(  X1 +   X2 +   X3 +   X4 +    X5-1 +     X6-0.5)+
#'     (W==2)*(2*X1 + 3*X2 +   X3 + 2*X4 + 2*(X5-1) + 2*(X6-0.5))+
#'     (W==3)*(3*X1 +   X2 + 2*X3 -   X4 -   (X5-1) -   (X6-0.5))+u
#'
#'   match1<-multilevelMatchX(Y,W,X)
#'   match2<-multilevelGPSMatch(Y,W,X,Trimming=FALSE,GPSM="multinomiallogisticReg")
#'   match3<-multilevelGPSMatch(Y,W,X,Trimming=TRUE,GPSM="multinomiallogisticReg")
#'   match4<-multilevelGPSStratification(Y,W,X,NS=10,GPSM="multinomiallogisticReg",linearp=0,nboot=50)
#'
#'   c(match1$tauestimate,match1$varestimate)
#'   c(match2$tauestimate,match2$varestimate)
#'   c(match3$tauestimate,match3$varestimate)
#'   c(match4$tauestimate,match4$varestimate)
#'
#' @import Matching boot nnet optmatch MASS
#'
#' @export
multilevelGPSStratification<-function(Y,W,X,NS,GPSM="multinomiallogisticReg",linearp=0,nboot){
  ### Generalized propensity score stratification
  ## PF.out=Propensity function estimation
  ## linearp=1 if use linear predictor for stratification
  ## X=design matrix
  ## Y=outcome
  ## W=treatment indicator
  ## nboot: the number of bootstrap replicates in variance estimator
  ## return values: point estimator,bootstrapping variance estimator
  #
  ## order the treatment increasingly

  N=length(Y) # number of observations
  X<-as.matrix(X)

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

  return(list(tauestimate=tauestimate,varestimate=bootvar))
}
