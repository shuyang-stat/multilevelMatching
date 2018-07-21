#' Matching on GPS with multilevel treatments
#'
#' @inheritParams multilevelMatchX
#' @param X A covariate matrix (p x n) with no intercept. When
#'   \code{GPSM="existing"}, then \code{X} must be a vector (1 x n) of
#'   user-specified propensity scores.
#' @param Trimming An indicator of whether trimming the sample to ensure overlap
#' @param GPSM An indicator of the methods used for estimating GPS, options
#'   include \code{"multinomiallogisticReg"}, \code{"ordinallogisticReg"} for
#'   proportional odds or cumulative logit model, and \code{"existing"} for
#'   user-specified propensity score via the parameter \code{X}. Defaults to
#'   \code{"multinomiallogisticReg"}
#'
#' @return A list element including:
#'   \itemize{
#'
#'     \item \code{tauestimate}: A vector of estimates for pairwise treatment
#'     effects
#'
#'     \item \code{varestimate}: A vector of variance estimates for
#'     \code{tauestimate}, using Abadie & Imbens (2006)'s method
#'
#'     \item \code{varestimateAI2012}: A vector of variance estimates for
#'     \code{tauestimate}, when matching on the generalized propensity score,
#'     using Abadie & Imbens (2016)'s method. This variance estimate takes into account
#'     of the uncertainty in estimating the GPS. This variable is named AI2012
#'     (not AI2016) for backwards compatibility.
#'
#'     \item \code{analysis_idx}: a list containing the indices_kept (analyzed)
#'     and indices_dropped (trimmed) based on Crump et al. (2009)'s method.
#'
#'   }
#'
#' @seealso \code{\link{multilevelMatchX}};
#'   \code{\link{multilevelGPSStratification}}
#'
#' @examples
#'   X <- c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9)
#'   Y <- c(102,105,120,130,100,80,94,108,96)
#'   W <- c(1,1,1,3,2,3,2,1,2)
#'   multilevelGPSMatch(Y,W,X,Trimming=0,GPSM="multinomiallogisticReg")
#'   multilevelGPSMatch(Y,W,X,Trimming=1,GPSM="multinomiallogisticReg")
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
#'   Crump, R. K., Hotz, V. J., Imbens, G. W., & Mitnik, O. A. (2009). Dealing
#'   with limited overlap in estimation of average treatment effects.
#'   Biometrika, 96(1), 187-199. \url{https://doi.org/10.1093/biomet/asn055}
#'
#' @export
multilevelGPSMatch <- function(
  Y,W,X,
  Trimming,
  GPSM = "multinomiallogisticReg"
){

  if ( Trimming == 1 ) {
    #PF modeling
    W.ref <- stats::relevel(as.factor(W),ref=1)
    temp <- utils::capture.output(PF.out <- nnet::multinom(W.ref~X))
    PF.fit <- stats::fitted(PF.out)
    ## identify sufficient overlap from Crump et al. (2009)'s method
    overlap.idx <- overlap(PF.fit)$idx

    W <- W[overlap.idx]
    X <- as.matrix(X)
    X <- X[overlap.idx,]
    Y <- Y[overlap.idx]
    analysisidx <- overlap.idx
  }
  if( Trimming == 0 ){
    analysisidx <- 1:length(Y)
  }

  ## order the treatment increasingly
  if ( 1-is.unsorted(W) ) {
    temp <- sort(W,index.return=TRUE)
    temp <- list(x=temp)
    temp$ix <- 1:length(W)
  }
  if(is.unsorted(W)){
    temp <- sort(W,index.return=TRUE)
  }
  W <- W[temp$ix]
  N <- length(Y) # number of observations
  X <- as.matrix(X)
  X <- X[temp$ix,]
  Y <- Y[temp$ix]


  trtnumber <- length(unique(W)) # number of treatment levels
  trtlevels <- unique(W) # all treatment levels
  pertrtlevelnumber <- table(W) # number of observations by treatment level
  taunumber <- trtnumber*(trtnumber+1)/2-trtnumber  # number of pairwise treatment effects


  # Modeling the generalized propensity score/function (PF)
  if (GPSM == "multinomiallogisticReg") {
    # message("Multinomial model fitted with reference level '1'")
    W.ref <- stats::relevel(as.factor(W), ref=1)
    temp <- utils::capture.output(PF.out <- nnet::multinom(W.ref~X))
    PF.fit <- stats::fitted(PF.out)
    vcov_coeff <- stats::vcov(PF.out)
  }
  if (GPSM == "ordinallogisticReg") {
    PF.out <- MASS::polr(as.factor(W)~X)
    PF.fit <- stats::fitted(PF.out)
  }
  if (GPSM == "existing") {
    # PF.out <- NULL
    ## Future debugging: defensive check that row sum of X equals 1
    PF.fit <- X
  }


  tauestimate <- varestimate <- varestimateAI2012 <- rep(NA,taunumber)
  meanw <- rep(NA,trtnumber)
  Yiw <- matrix(NA,N,trtnumber) #Yiw is the full imputed data set
  Kiw <- sigsqiw <- matrix(NA,N,1)  #Kiw is vector of number of times unit i used as a match

  Matchmat <- matrix(NA,N,trtnumber*2)
  cname <- c()
  for ( kk in 1:trtnumber ) {
    thistrt <- trtlevels[kk]
    cname <- c(cname,c(paste(paste(paste("m",thistrt,sep=""),".",sep=""),1,sep=""),
                     paste(paste(paste("m",thistrt,sep=""),".",sep=""),2,sep="")))
  }
  colnames(Matchmat) <- cname


  for ( kk in 1:trtnumber ) {
    thistrt <- trtlevels[kk]
    if ( kk == 1 ) {fromto <- 1:pertrtlevelnumber[1]}
    if ( kk > 1 ) {fromto <- (1:pertrtlevelnumber[kk])+sum(pertrtlevelnumber[1:(kk-1)])}
    W1 <- W!=thistrt
    out1 <- Matching::Match(Y=Y,Tr=W1,X=PF.fit[,kk],distance.tolerance=0,ties=FALSE,Weight=2)
    mdata1 <- out1$mdata
    meanw[kk] <- stats::weighted.mean(c(Y[which(W==thistrt)],mdata1$Y[which(mdata1$Tr==0)]),c(rep(1,length(which(W==thistrt))),out1$weights))
    Kiw[fromto,1] <- table(factor(out1$index.control,levels=fromto))
    Yiw[which(W==thistrt),kk] <- Y[which(W==thistrt)]
    Yiw[which(W!=thistrt),kk] <- mdata1$Y[which(mdata1$Tr==0)]

    WW1 <- W==thistrt
    out11 <- Matching::Match(Y=rep(Y[which(WW1)],times=2),Tr=rep(c(1,0),each=sum(WW1)),
                 X=c(PF.fit[which(WW1),kk],PF.fit[which(WW1),kk]),M=1,distance.tolerance=0,ties=FALSE,Weight=2,
                 restrict=matrix(c(1:sum(WW1),(1:sum(WW1))+sum(WW1),rep(-1,sum(WW1))),nrow=sum(WW1),ncol=3,byrow=FALSE))

    mdata11 <- out11$mdata
    temp11 <- (mdata11$Y[which(mdata11$Tr==1)]-mdata11$Y[which(mdata11$Tr==0)])^2/2
    sigsqiw[which(W==thistrt),1] <- temp11

    thiscnames <- c(paste(paste(paste("m",thistrt,sep=""),".",sep=""),1,sep=""),
                  paste(paste(paste("m",thistrt,sep=""),".",sep=""),2,sep=""))

    # find two outsiders closest
    findmatch1 <- Matching::Match(Y=Y,Tr=W1,X=PF.fit[,kk],distance.tolerance=0,ties=FALSE,Weight=2,M=2)
    Matchmat[unique(findmatch1$index.treated),thiscnames] <- matrix(findmatch1$index.control,ncol=2,byrow=TRUE)
    # find one insider closest
    out111 <- Matching::Match(Y=rep(Y[which(WW1)],times=2),Tr=rep(c(0,1),each=sum(WW1)),
                  X=c(PF.fit[which(WW1),kk],PF.fit[which(WW1),kk]),M=1,distance.tolerance=0,ties=FALSE,Weight=2,
                  restrict=matrix(c(1:sum(WW1),(1:sum(WW1))+sum(WW1),rep(-1,sum(WW1))),nrow=sum(WW1),ncol=3,byrow=FALSE))
    Matchmat[which(WW1),thiscnames] <- matrix(c(which(WW1),which(WW1)[out111$index.control]),ncol=2,byrow=FALSE)

  }

  cnt <- 0
  cname1 <- c()
  for ( jj in 1:(trtnumber-1) ) {
    for ( kk in (jj+1):trtnumber ) {
      cnt <- cnt+1
      thistrt <- trtlevels[jj]
      thattrt <- trtlevels[kk]
      cname1 <- c(cname1,paste(paste(paste(paste(paste("EY(",thattrt,sep=""),")",sep=""),"-EY(",sep=""),thistrt,sep=""),")",sep=""))
      tauestimate[cnt] <- meanw[kk]-meanw[jj]
      varestimate[cnt] <- mean((Yiw[,kk]-Yiw[,jj]-(meanw[kk]-meanw[jj]))^2)+mean((Kiw^2+Kiw)*sigsqiw*(W==thistrt | W==thattrt))
    }
  }
  varestimate <- varestimate/N
  names(tauestimate) <- cname1
  names(varestimate) <- cname1
  names(varestimateAI2012) <- cname1

  if (GPSM=="multinomiallogisticReg") {

    I.inv <- vcov_coeff
    ## Adjustment term c'(I^-1)c
    X <- as.matrix(X)
    Cmat <- matrix(0,N,(dim(X)[2]+1)*(trtnumber-1))
    Cvec <- matrix(0,trtnumber,(dim(X)[2]+1)*(trtnumber-1))
    for(kkk in 1:trtnumber){
      thistrt <- trtlevels[kkk]
      thiscnames <- c(paste(paste(paste("m",thistrt,sep=""),".",sep=""),1,sep=""),
                    paste(paste(paste("m",thistrt,sep=""),".",sep=""),2,sep=""))
      Y11 <- matrix(Y[Matchmat[,c(thiscnames)]],ncol=2,byrow=FALSE)
      mY11 <- apply(Y11,1,mean)
      for (kk in 1:(trtnumber-1) ) {
        for (jj in 1:(dim(X)[2]+1) ) {
          if (jj == 1) {}
          if (jj > 1) {
            X11 <- matrix(X[Matchmat[,c(thiscnames)],(jj-1)],ncol=2,byrow=FALSE)
            mX11 <- apply(X11,1,mean)
            C1.X1Y <- apply((X11-mX11)*(Y11-mY11),1,sum)
            if(kkk==(kk+1)){C1.X1Y <- C1.X1Y*(1-PF.fit[,kk+1])}
            else if(kkk!=(kk+1))C1.X1Y <- C1.X1Y*(-PF.fit[,kk+1])
            Cmat[,(dim(X)[2]+1)*(kk-1)+jj] <- C1.X1Y
          }
        }
      }
      Cvec[kkk,] <- apply(Cmat,2,mean)
    }

    for ( jj in 1:(trtnumber-1) ) {
      for ( kk in (jj+1):trtnumber ) {
        thistrt <- trtlevels[jj]
        thattrt <- trtlevels[kk]
        cname1 <- c(paste(paste(paste(paste(paste("EY(",thattrt,sep=""),")",sep=""),"-EY(",sep=""),thistrt,sep=""),")",sep=""))
        varestimateAI2012[cname1] <- varestimate[cname1]-
          t(Cvec[jj,]+Cvec[kk,])%*%vcov_coeff%*%(Cvec[jj,]+Cvec[kk,])
      }
    }
  }

  out_list <- list(
    tauestimate = tauestimate,
    varestimate = varestimate,
    varestimateAI2012 = varestimateAI2012,
    analysisidx = analysisidx
  )

  return(out_list)
}
