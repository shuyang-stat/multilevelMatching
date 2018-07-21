
# #' Statistic for bootstrapping the standard errors from GPS Stratification
# #'
# #'
estforboot <- function(
  data,indices,GPSM,linearp,trtnumber,trtlevels,taunumber,NS
){
  ###boot function
  #
  d <- data[indices,] # allows boot to select sample
  W<-d[,"W"]
  W<-as.factor(W)
  Y<-d[,"Y"]
  X<-d[,-c(1:2)]

  #PF modeling
  if(GPSM=="multinomiallogisticReg"){
    W.ref <- stats::relevel(as.factor(W),ref="1")
    temp <- utils::capture.output(PF.out <- nnet::multinom(W.ref~X))
    PF.fit <- stats::fitted(PF.out)
    if(linearp==1){
      beta <- stats::coef(PF.out)
      Xbeta<-X%*%t(beta[,-1])
      PF.fit[,-1]<-Xbeta
    }
  }
  if(GPSM=="ordinallogisticReg"){
    PF.out <- MASS::polr(as.factor(W)~X)
    PF.fit <- stats::fitted(PF.out)
  }
  if(GPSM=="existing"){
    PF.fit <- X
  }

  meanwj<-numberwj<-matrix(NA,trtnumber,NS)

  for(kk in 1:trtnumber){
    pwx<-PF.fit[,kk]
    ranking <- stats::ave( pwx, FUN=function(x) {
      cut(x,stats::quantile(pwx,(0:NS)/NS,type=2),
          include.lowest=TRUE,right=FALSE)
    })
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
