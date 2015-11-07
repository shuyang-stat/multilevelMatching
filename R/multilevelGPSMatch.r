#' Matching on GPS with multilevel treatments
#'
#' @param Y a continuous response vector (1 x n)
#' @param W a treatment vector (1 x n) with numerical values indicating treatment groups
#' @param X a covariate matrix (p x n) with no intercept
#' @param Trimming an indicator of whether trimming the sample to ensure overlap
#' @param GPSM an indicator of the methods used for estimating GPS, options include "multinomiallogisticReg", "ordinallogisticReg", and "existing"
#'
#' @return
#' \itemize{
#'
#'  \item tauestimate:  a vector of estimates for pairwise treatment effects
#'
#'  \item varestimate:  a vector of variance estimates for tauestimate, using Abadie&Imbens(2006)'s method
#'
#'  \item varestimateAI2012:  a vector of variance estimates for tauestimate, when matching on the generalized propensity score, using Abadie&Imbens(2012)'s method.
#'  This variance estimate takes into account of the uncertainty in estimating the GPS.
#'
#'  \item analysisidx: the index of units after trimming based on Crump et al. (2009)'s method.
#'
#' }
#'
#' @seealso \code{\link{multilevelMatchX}}; \code{\link{multilevelGPSStratification}}
#'
#' @examples
#'   X<-c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9)
#'   Y<-c(102,105,120,130,100,80,94,108,96)
#'   W<-c(1,1,1,3,2,3,2,1,2)
#'   multilevelGPSMatch(Y,W,X,Trimming=0,GPSM="multinomiallogisticReg")
#'   multilevelGPSMatch(Y,W,X,Trimming=1,GPSM="multinomiallogisticReg")
#'
#' @import Matching boot nnet optmatch MASS
#'
#' @export
multilevelGPSMatch<-function(Y,W,X,Trimming,GPSM="multinomiallogisticReg"){

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
