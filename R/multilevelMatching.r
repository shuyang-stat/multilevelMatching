#' multilevel treatment matching on X
#'
#' @param Y continuous response vector (1 x n)
#' @param W treatment vector (1 x n) with numerical values indicating treatment groups
#' @param X covariate matrix (p x n) with no intercept
#' @param Trimming
#' @param GPSM indicator of the methods used for estimating GPS, options include "multinomiallogisticReg", "ordinallogisticReg", and "existing"
#' @param NS (only required in the function multilevelGPSStratification) number of strata
#' @param linearp (only required in the function multilevelGPSStratification) indicator of subclassification on GPS (=0) or linear predictor of GPS (=1)
#' @param nboot (only required in the function multilevelGPSStratification) number of boot replicates for variance estimation
#'
#'
#' @return A list with 2 elements: tauestimate, varestimate, varestimateAI2012(only for multilevelGPSMatch), analysisidx (when Trimming=1)
#' tauestimate is a vector of estimates for pairwise treatment effects
#' varestimate is a vector of variance estimates for tauestimate, using Abadie&Imbens(2006)'s method
#' varestimateAI2012 is a vector of variance estimates for tauestimate,
#'  when matching on the generalized propensity score, using Abadie&Imbens(2012)'s method.
#'  This variance estimate takes into account of the uncertainty in estimating the GPS.
#' analysisidx is the index of units after trimming based on Crump et al. (2009)'s method.
#'
#' @examples
#'   X<-c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9)
#'   Y<-c(102,105,120,130,100,80,94,108,96)
#'   W<-c(1,1,1,3,2,3,2,1,2)
#'   multilevelMatchX(Y,W,X)
#'   multilevelGPSMatch(Y,W,X,Trimming=0,GPSM="multinomiallogisticReg")
#'   multilevelGPSMatch(Y,W,X,Trimming=1,GPSM="multinomiallogisticReg")
#'
#'  set.seed(111)
#'   n    <- 10000*6
#'   # X1-X3 3 MVN var 2, 1, 1, covars 1, -1, -.5
#'   vars   <- c(2,1,1)
#'   covars <- c(1,-1,-.5)
#'   mu     <- c(0,0,0)
#'   tau    <- 1
#'   Sigma <- diag(vars)
#'   Sigma[2,1] <- Sigma[1,2] <- covars[1]
#'   Sigma[3,1] <- Sigma[1,3] <- covars[2]
#'   Sigma[3,2] <- Sigma[2,3] <- covars[3]
#'   trt1 <- 500; trt1
#'   trt2 <- 500; trt2
#'   trt3 <- 500; trt3
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
#'   TX<-X
#'   mattau2<-mean (2*TX[,1] + 3*TX[,2] +   TX[,3] + 2*TX[,4] + 2*(TX[,5]-1) + 2*(TX[,6]-0.5))-mean(  TX[,1] +   TX[,2] + TX[,3] +   TX[,4] +    TX[,5]-1 +     TX[,6]-0.5)
#'   mattau3<-mean (3*TX[,1] +   TX[,2] + 2*TX[,3] -   TX[,4] -   (TX[,5]-1) -   (TX[,6]-0.5))-mean(  TX[,1] +   TX[,2] + TX[,3] +   TX[,4] +    TX[,5]-1 +     TX[,6]-0.5)
#'   mattau23<-mean(3*TX[,1] +   TX[,2] + 2*TX[,3] -   TX[,4] -   (TX[,5]-1) -   (TX[,6]-0.5))-mean(2*TX[,1] + 3*TX[,2] + TX[,3] + 2*TX[,4] + 2*(TX[,5]-1) + 2*(TX[,6]-0.5))
#'
#'   idp<-which(W!=3)
#'   TX<-X[idp,]
#'   pair2<-mean(2*TX[,1] + 3*TX[,2] + TX[,3] + 2*TX[,4] + 2*(TX[,5]-1) + 2*(TX[,6]-0.5))-mean(TX[,1] + TX[,2] + TX[,3] + TX[,4] + TX[,5]-1 + TX[,6]-0.5)
#'   idp<-which(W!=2)
#'   TX<-X[idp,]
#'   pair3<-mean(3*TX[,1] + TX[,2] + 2*TX[,3] - TX[,4] - (TX[,5]-1) - (TX[,6]-0.5))-mean(TX[,1] + TX[,2] + TX[,3] + TX[,4] + TX[,5]-1 + TX[,6]-0.5)
#'   idp<-which(W!=1)
#'   TX<-X[idp,]
#'   pair23<-mean(3*TX[,1] + TX[,2] + 2*TX[,3] - TX[,4] - (TX[,5]-1) - (TX[,6]-0.5))-mean(2*TX[,1] + 3*TX[,2] + TX[,3] + 2*TX[,4] + 2*(TX[,5]-1) + 2*(TX[,6]-0.5))
#'
#'   estimand<-c(mattau2,mattau3,mattau23)
#'   match1<-multilevelMatchX(Y,W,X)
#'   match2<-multilevelGPSMatch(Y,W,X,Trimming=FALSE,GPSM="multinomiallogisticReg")
#'   match3<-multilevelGPSMatch(Y,W,X,Trimming=TRUE,GPSM="multinomiallogisticReg")
#'   match4<-multilevelGPSStratification(Y,W,X,NS=10,GPSM="multinomiallogisticReg",linearp=0,nboot=50)
#'
#'   estimand
#'   c(match1$tauestimate,match1$varestimate)
#'   c(match2$tauestimate,match2$varestimate)
#'   c(match3$tauestimate,match3$varestimate)
#'   c(match4$tauestimate,match4$varestimate)
#'
#' @export multilevelMatchX multilevelGPSMatch  multilevelGPSStratification
#' @seealso \code{\link{multilevelMatchX}},
#' \code{\link{multilevelGPSMatch}},
#' \code{\link{multilevelGPSStratification}}
#' \code{\link{overlap}}
#' \code{\link{estforboot}}
#'
multilevelMatchX<-function(Y,W,X){

  ## order the treatment increasingly
  if(1-is.unsorted(W)){
    temp<-sort(W,index.return=TRUE)
    temp<-list(x=temp)
    temp$ix<-1:length(W)
  }
  if(is.unsorted(W)){
    temp<-sort(W,index.return=TRUE)
  }
  W<-W[temp$ix]
  N=length(Y) # number of observations

  X<-as.matrix(X)
  X<-X[temp$ix,]
  Y<-Y[temp$ix]

  trtnumber<-length(unique(W)) # number of treatment levels
  trtlevels<-unique(W) # all treatment levels
  pertrtlevelnumber<-table(W) # number of observations by treatment level
  taunumber<-trtnumber*(trtnumber+1)/2-trtnumber  # number of pairwise treatment effects

  tauestimate<-varestimate<-rep(NA,taunumber)
  meanw<-rep(NA,trtnumber)

  Yiw<-matrix(NA,N,trtnumber) #Yiw is the full imputed data set
  Kiw<-sigsqiw<-matrix(NA,N,1)  #Kiw is vector of number of times unit i used as a match

  Matchmat<-matrix(NA,N,trtnumber*2)
  cname<-c()
  for(kk in 1:trtnumber){
    thistrt<-trtlevels[kk]
    cname<-c(cname,c(paste(paste(paste("m",thistrt,sep=""),".",sep=""),1,sep=""),
                     paste(paste(paste("m",thistrt,sep=""),".",sep=""),2,sep="")))
  }
  colnames(Matchmat)<-cname


  for(kk in 1:trtnumber){

    thistrt<-trtlevels[kk]
    if(kk==1){fromto<-1:pertrtlevelnumber[1]}
    if(kk>1){fromto<-(1:pertrtlevelnumber[kk])+sum(pertrtlevelnumber[1:(kk-1)])}
    W1<-W!=thistrt
    out1 <- Match(Y=Y,Tr=W1,X=X,distance.tolerance=0,ties=FALSE,Weight=2)
    mdata1<-out1$mdata
    meanw[kk]<-weighted.mean(c(Y[which(W==thistrt)],mdata1$Y[which(mdata1$Tr==0)]),c(rep(1,length(which(W==thistrt))),out1$weights))
    Kiw[fromto,1]<-table(factor(out1$index.control,levels=fromto))
    Yiw[which(W==thistrt),kk]<- Y[which(W==thistrt)]
    Yiw[which(W!=thistrt),kk]<-mdata1$Y[which(mdata1$Tr==0)]
    WW1<-W==thistrt
    X<-as.matrix(X)
    out11<-Match(Y=rep(Y[which(WW1)],times=2),Tr=rep(c(1,0),each=sum(WW1)),
                 X=rbind(as.matrix(X[which(WW1),]),as.matrix(X[which(WW1),])),M=1,distance.tolerance=0,ties=FALSE,Weight=2,
                 restrict=matrix(c(1:sum(WW1),(1:sum(WW1))+sum(WW1),rep(-1,sum(WW1))),nrow=sum(WW1),ncol=3,byrow=FALSE))

    mdata11<-out11$mdata
    temp11<-(mdata11$Y[which(mdata11$Tr==1)]-mdata11$Y[which(mdata11$Tr==0)])^2/2
    sigsqiw[which(W==thistrt),1]<-temp11

    thiscnames<-c(paste(paste(paste("m",thistrt,sep=""),".",sep=""),1,sep=""),
                  paste(paste(paste("m",thistrt,sep=""),".",sep=""),2,sep=""))

    # find two outsiders closest
    findmatch1<-Match(Y=Y,Tr=W1,X=X,distance.tolerance=0,ties=FALSE,Weight=2,M=2)
    Matchmat[unique(findmatch1$index.treated),thiscnames]<-matrix(findmatch1$index.control,ncol=2,byrow=TRUE)
    # find one insider closest
    out111<-Match(Y=rep(Y[which(WW1)],times=2),Tr=rep(c(0,1),each=sum(WW1)),
                  X=rbind(as.matrix(X[which(WW1),]),as.matrix(X[which(WW1),])),M=1,distance.tolerance=0,ties=FALSE,Weight=2,
                  restrict=matrix(c(1:sum(WW1),(1:sum(WW1))+sum(WW1),rep(-1,sum(WW1))),nrow=sum(WW1),ncol=3,byrow=FALSE))
    Matchmat[which(WW1),thiscnames]<-matrix(c(which(WW1),which(WW1)[out111$index.control]),ncol=2,byrow=FALSE)

  }
  cnt<-0
  cname1<-c()
  for(jj in 1:(trtnumber-1)){
    for(kk in (jj+1):trtnumber){
      cnt<-cnt+1
      thistrt<-trtlevels[jj]
      thattrt<-trtlevels[kk]
      cname1<-c(cname1,paste(paste(paste(paste(paste("EY(",thattrt,sep=""),")",sep=""),"-EY(",sep=""),thistrt,sep=""),")",sep=""))
      tauestimate[cnt]<-meanw[kk]-meanw[jj]
      varestimate[cnt]<-mean((Yiw[,kk]-Yiw[,jj]-(meanw[kk]-meanw[jj]))^2)+mean((Kiw^2+Kiw)*sigsqiw*(W==thistrt | W==thattrt))
    }
  }
  varestimate<-varestimate/N
  names(tauestimate)<-cname1
  names(varestimate)<-cname1

  return(list(tauestimate=tauestimate,varestimate=varestimate))

}


multilevelGPSMatch<-function(Y,W,X,Trimming,GPSM="multinomiallogisticReg"){
  ###GPSM="multinomiallogisticReg"/"ordinallogisticReg"/"existing"(X is the fitted GPS by users)
  ##return values: point estimator; A&I(2006)'s variance estimator; A&I(2012)'s variance estimator
  #
  ###### Trimming ##################################################

  if(Trimming==1){
    #PF modeling
    W.ref <- relevel(as.factor(W),ref=1)
    temp<-capture.output(PF.out <- multinom(W.ref~X))
    PF.fit <- fitted(PF.out)
    ## identify sufficient overlap
    overlap.idx<-overlap(PF.fit)$idx

    W <- W[overlap.idx]
    X <- as.matrix(X)
    X <- X[overlap.idx,]
    Y <- Y[overlap.idx]
    analysisidx<-overlap.idx
  }
  if(Trimming==0){
    analysisidx<-1:length(Y)
  }
  ## order the treatment increasingly
  if(1-is.unsorted(W)){
    temp<-sort(W,index.return=TRUE)
    temp<-list(x=temp)
    temp$ix<-1:length(W)
  }
  if(is.unsorted(W)){
    temp<-sort(W,index.return=TRUE)
  }
  W<-W[temp$ix]
  N=length(Y) # number of observations
  X<-as.matrix(X)
  X<-X[temp$ix,]
  Y<-Y[temp$ix]


  trtnumber<-length(unique(W)) # number of treatment levels
  trtlevels<-unique(W) # all treatment levels
  pertrtlevelnumber<-table(W) # number of observations by treatment level
  taunumber<-trtnumber*(trtnumber+1)/2-trtnumber  # number of pairwise treatment effects


  #PF modeling
  if(GPSM=="multinomiallogisticReg"){
    W.ref <- relevel(as.factor(W),ref=1)
    temp<-capture.output(PF.out <- multinom(W.ref~X))
    PF.fit <- fitted(PF.out)
    vcov_coeff <- vcov(PF.out)
  }
  if(GPSM=="ordinallogisticReg"){
    PF.out <- polr(as.factor(W)~X)
    PF.fit <- fitted(PF.out)
  }
  if(GPSM=="existing"){
    #need to check the row sum of X is 1 - debug1
    PF.fit <- X
  }

  tauestimate<-varestimate<-varestimateAI2012<-rep(NA,taunumber)
  meanw<-rep(NA,trtnumber)

  Yiw<-matrix(NA,N,trtnumber) #Yiw is the full imputed data set
  Kiw<-sigsqiw<-matrix(NA,N,1)  #Kiw is vector of number of times unit i used as a match

  Matchmat<-matrix(NA,N,trtnumber*2)
  cname<-c()
  for(kk in 1:trtnumber){
    thistrt<-trtlevels[kk]
    cname<-c(cname,c(paste(paste(paste("m",thistrt,sep=""),".",sep=""),1,sep=""),
                     paste(paste(paste("m",thistrt,sep=""),".",sep=""),2,sep="")))
  }
  colnames(Matchmat)<-cname


  for(kk in 1:trtnumber){
    thistrt<-trtlevels[kk]
    if(kk==1){fromto<-1:pertrtlevelnumber[1]}
    if(kk>1){fromto<-(1:pertrtlevelnumber[kk])+sum(pertrtlevelnumber[1:(kk-1)])}
    W1<-W!=thistrt
    out1 <- Match(Y=Y,Tr=W1,X=PF.fit[,kk],distance.tolerance=0,ties=FALSE,Weight=2)
    mdata1<-out1$mdata
    meanw[kk]<-weighted.mean(c(Y[which(W==thistrt)],mdata1$Y[which(mdata1$Tr==0)]),c(rep(1,length(which(W==thistrt))),out1$weights))
    Kiw[fromto,1]<-table(factor(out1$index.control,levels=fromto))
    Yiw[which(W==thistrt),kk]<- Y[which(W==thistrt)]
    Yiw[which(W!=thistrt),kk]<-mdata1$Y[which(mdata1$Tr==0)]

    WW1<-W==thistrt
    out11<-Match(Y=rep(Y[which(WW1)],times=2),Tr=rep(c(1,0),each=sum(WW1)),
                 X=c(PF.fit[which(WW1),kk],PF.fit[which(WW1),kk]),M=1,distance.tolerance=0,ties=FALSE,Weight=2,
                 restrict=matrix(c(1:sum(WW1),(1:sum(WW1))+sum(WW1),rep(-1,sum(WW1))),nrow=sum(WW1),ncol=3,byrow=FALSE))

    mdata11<-out11$mdata
    temp11<-(mdata11$Y[which(mdata11$Tr==1)]-mdata11$Y[which(mdata11$Tr==0)])^2/2
    sigsqiw[which(W==thistrt),1]<-temp11

    thiscnames<-c(paste(paste(paste("m",thistrt,sep=""),".",sep=""),1,sep=""),
                  paste(paste(paste("m",thistrt,sep=""),".",sep=""),2,sep=""))

    # find two outsiders closest
    findmatch1<-Match(Y=Y,Tr=W1,X=PF.fit[,kk],distance.tolerance=0,ties=FALSE,Weight=2,M=2)
    Matchmat[unique(findmatch1$index.treated),thiscnames]<-matrix(findmatch1$index.control,ncol=2,byrow=TRUE)
    # find one insider closest
    out111<-Match(Y=rep(Y[which(WW1)],times=2),Tr=rep(c(0,1),each=sum(WW1)),
                  X=c(PF.fit[which(WW1),kk],PF.fit[which(WW1),kk]),M=1,distance.tolerance=0,ties=FALSE,Weight=2,
                  restrict=matrix(c(1:sum(WW1),(1:sum(WW1))+sum(WW1),rep(-1,sum(WW1))),nrow=sum(WW1),ncol=3,byrow=FALSE))
    Matchmat[which(WW1),thiscnames]<-matrix(c(which(WW1),which(WW1)[out111$index.control]),ncol=2,byrow=FALSE)

  }

  cnt<-0
  cname1<-c()
  for(jj in 1:(trtnumber-1)){
    for(kk in (jj+1):trtnumber){
      cnt<-cnt+1
      thistrt<-trtlevels[jj]
      thattrt<-trtlevels[kk]
      cname1<-c(cname1,paste(paste(paste(paste(paste("EY(",thattrt,sep=""),")",sep=""),"-EY(",sep=""),thistrt,sep=""),")",sep=""))
      tauestimate[cnt]<-meanw[kk]-meanw[jj]
      varestimate[cnt]<-mean((Yiw[,kk]-Yiw[,jj]-(meanw[kk]-meanw[jj]))^2)+mean((Kiw^2+Kiw)*sigsqiw*(W==thistrt | W==thattrt))
    }
  }
  varestimate<-varestimate/N
  names(tauestimate)<-cname1
  names(varestimate)<-cname1
  names(varestimateAI2012)<-cname1


  if(GPSM=="multinomiallogisticReg"){
    I.inv<-vcov_coeff
    ## Adjustment term c'(I^-1)c
    X<-as.matrix(X)
    Cmat<-matrix(0,N,(dim(X)[2]+1)*(trtnumber-1))
    Cvec<-matrix(0,trtnumber,(dim(X)[2]+1)*(trtnumber-1))
    for(kkk in 1:trtnumber){
      thistrt<-trtlevels[kkk]
      thiscnames<-c(paste(paste(paste("m",thistrt,sep=""),".",sep=""),1,sep=""),
                    paste(paste(paste("m",thistrt,sep=""),".",sep=""),2,sep=""))
      Y11<-matrix(Y[Matchmat[,c(thiscnames)]],ncol=2,byrow=FALSE)
      mY11<-apply(Y11,1,mean)
      for(kk in 1:(trtnumber-1)){
        for(jj in 1:(dim(X)[2]+1)){
          if(jj==1){}
          if(jj>1){
            X11<-matrix(X[Matchmat[,c(thiscnames)],(jj-1)],ncol=2,byrow=FALSE)
            mX11<-apply(X11,1,mean)
            C1.X1Y<-apply((X11-mX11)*(Y11-mY11),1,sum)
            C1.X1Y<-C1.X1Y*(-PF.fit[,kk+1])
            Cmat[,(dim(X)[2]+1)*(kk-1)+jj]<-C1.X1Y
          }
        }
      }
      Cvec[kkk,]<-apply(Cmat,2,mean)
    }

    for(jj in 1:(trtnumber-1)){
      for(kk in (jj+1):trtnumber){
        thistrt<-trtlevels[jj]
        thattrt<-trtlevels[kk]
        cname1<-c(cname1,paste(paste(paste(paste(paste("EY(",thattrt,sep=""),")",sep=""),"-EY(",sep=""),thistrt,sep=""),")",sep=""))
        varestimateAI2012[cname1]<-varestimate[cname1]-
          t(Cvec[jj,]+Cvec[kk,])%*%vcov_coeff%*%(Cvec[jj,]+Cvec[kk,])
      }
    }
  }
  return(list(tauestimate=tauestimate,
              varestimate=varestimate,
              varestimateAI2012=varestimateAI2012,
              analysisidx=analysisidx))
}


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

estforboot <- function(data,indices,GPSM,linearp,trtnumber,trtlevels,taunumber,NS) {
  ###boot function
  #
  d <- data[indices,] # allows boot to select sample
  W<-d[,"W"]
  W<-as.factor(W)
  Y<-d[,"Y"]
  X<-d[,-c(1:2)]

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
      numberwj[kk,jj]<-sum(ranking==jj)
    }
    numberwj<-numberwj*(1-is.na(meanwj))
    meanw<-apply(meanwj*numberwj,1,sum,na.rm=TRUE)/apply(numberwj,1,sum)
  }

  tauestimate<-rep(NA,taunumber)
  cnt<-0
  for(jj in 1:(trtnumber-1)){
    for(kk in (jj+1):trtnumber){
      cnt<-cnt+1
      tauestimate[cnt]<-meanw[kk]-meanw[jj]
    }
  }

  return(tauestimate)
}

overlap<-function(PF.fit){
  obj<-function(alpha){
    gx<-apply(1/PF.fit,1,sum)
    id<-(gx<alpha)
    mean(gx*id)/(mean(id)^2)
  }
  gx<-apply(1/PF.fit,1,sum)
  #alpha1<-optim(median(gx), obj)$par
  suppressWarnings( alpha1<-optim(median(gx), obj)$par)
  alpha<-2*obj(alpha1)
  idx<-which(gx<alpha)
  return(list(alpha=alpha,idx=idx))
}






