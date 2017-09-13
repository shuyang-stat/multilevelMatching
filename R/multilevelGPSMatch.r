#' Matching on GPS with multilevel treatments
#'
#' @param Y a continuous response vector (1 x n)
#' @param W a treatment vector (1 x n) with numerical values indicating
#'   treatment groups
#' @param X a covariate matrix (p x n) with no intercept. When GPSM="existing",
#'   then X must be a vector (1 x n) of user-specified propensity scores.
#' @param Trimming an indicator of whether trimming the sample to ensure overlap
#' @param GPSM an indicator of the methods used for estimating GPS, options
#'   include "multinomiallogisticReg", "ordinallogisticReg" for proportional
#'   odds or cumulative logit model, and "existing" for user-specified
#'   propensity score via the parameter X.
#' @param model_options A list of the options to pass to propensity model.
#'   Currently under development. Can only pass reference level to multinomial
#'   logisitc regression.
#'
#' @return according to \code{\link{estimateTau}}, including at most: \itemize{
#'
#'   \item tauestimate:  a vector of estimates for pairwise treatment effects
#'
#'   \item varestimate:  a vector of variance estimates for tauestimate, using
#'   Abadie&Imbens(2006)'s method
#'
#'   \item varestimateAI2012:  a vector of variance estimates for tauestimate,
#'   when matching on the generalized propensity score, using
#'   Abadie&Imbens(2012)'s method. This variance estimate takes into account of
#'   the uncertainty in estimating the GPS.
#'
#'   \item analysis_idx: a list containing the indices_kept (analyzed) and
#'   indices_dropped (trimmed) based on Crump et al. (2009)'s method.
#'
#'   }
#'
#' @seealso \code{\link{multilevelMatchX}};
#'   \code{\link{multilevelGPSStratification}}
#'
#' @examples
#'   X<-c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9)
#'   Y<-c(102,105,120,130,100,80,94,108,96)
#'   W<-c(1,1,1,3,2,3,2,1,2)
#'   multilevelGPSMatch(Y,W,X,Trimming=0,GPSM="multinomiallogisticReg")
#'   multilevelGPSMatch(Y,W,X,Trimming=1,GPSM="multinomiallogisticReg")
#'
#' @import Matching boot nnet MASS
#'
#' @export
multilevelGPSMatch <- function(Y,W,X,Trimming,GPSM="multinomiallogisticReg",
                               model_options = list(reference_level = "1")){

  match_method <- "MatchOnGPS"

  #probably move this into prepareData
  if (!is.null(model_options)){
    if (!is.list(model_options)){
      stop("model_options must be a list or NULL")
    }
    # if (GPSM!= "existing") {
    #
    # }
    if (GPSM == "multinomiallogisticReg"){
      if (!"reference_level" %in% names(model_options)){
        stop("User must supply model_options$reference_level")
      } ## defensive programming for correct level variable?
    }
  }

  prepared_data <- prepareData(
    Y=Y, W=W, X=X,
    match_method = match_method,
    Trimming = Trimming
    #Trimming_fit_args
  )
  W <- prepared_data$W
  X <- prepared_data$X
  Y <- prepared_data$Y
  N <- prepared_data$N
  trtnumber <- prepared_data$trtnumber
  trtlevels <- prepared_data$trtlevels
  pertrtlevelnumber <- prepared_data$pertrtlevelnumber
  taunumber <- prepared_data$taunumber
  analysis_idx <- prepared_data$analysis_idx

  #PF modeling
  if (GPSM == "multinomiallogisticReg") {
    # message("Multinomial model fitted with reference level '1'")
    ##This should probably be user-specified via the dots argument
    W.ref <- stats::relevel(as.factor(W),ref=model_options$reference_level)
    temp <- utils::capture.output(PF.out <- nnet::multinom(W.ref~X))
    PF.fit <- stats::fitted(PF.out)
    vcov_coeff <- stats::vcov(PF.out)
  }
  if (GPSM == "ordinallogisticReg") {
    PF.out <- MASS::polr(as.factor(W)~X)
    PF.fit <- stats::fitted(PF.out)
  }
  if (GPSM == "existing") {
    ## bug checks migrated to prepareData()
    PF.fit <- X
  }


  # tauestimate<-varestimate<-varestimateAI2012<-rep(NA,taunumber)
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
    out1 <- Matching::Match(Y=Y,Tr=W1,X=PF.fit[,kk],distance.tolerance=0,ties=FALSE,Weight=2)
    mdata1<-out1$mdata
    meanw[kk]<-stats::weighted.mean(c(Y[which(W==thistrt)],mdata1$Y[which(mdata1$Tr==0)]),c(rep(1,length(which(W==thistrt))),out1$weights))
    Kiw[fromto,1]<-table(factor(out1$index.control,levels=fromto))
    Yiw[which(W==thistrt),kk]<- Y[which(W==thistrt)]
    Yiw[which(W!=thistrt),kk]<-mdata1$Y[which(mdata1$Tr==0)]

    WW1<-W==thistrt
    out11<-Matching::Match(Y=rep(Y[which(WW1)],times=2),Tr=rep(c(1,0),each=sum(WW1)),
                 X=c(PF.fit[which(WW1),kk],PF.fit[which(WW1),kk]),M=1,distance.tolerance=0,ties=FALSE,Weight=2,
                 restrict=matrix(c(1:sum(WW1),(1:sum(WW1))+sum(WW1),rep(-1,sum(WW1))),nrow=sum(WW1),ncol=3,byrow=FALSE))

    mdata11<-out11$mdata
    temp11<-(mdata11$Y[which(mdata11$Tr==1)]-mdata11$Y[which(mdata11$Tr==0)])^2/2
    sigsqiw[which(W==thistrt),1]<-temp11

    thiscnames<-c(paste(paste(paste("m",thistrt,sep=""),".",sep=""),1,sep=""),
                  paste(paste(paste("m",thistrt,sep=""),".",sep=""),2,sep=""))

    # find two outsiders closest
    findmatch1<-Matching::Match(Y=Y,Tr=W1,X=PF.fit[,kk],distance.tolerance=0,ties=FALSE,Weight=2,M=2)
    Matchmat[unique(findmatch1$index.treated),thiscnames]<-matrix(findmatch1$index.control,ncol=2,byrow=TRUE)
    # find one insider closest
    out111<-Matching::Match(Y=rep(Y[which(WW1)],times=2),Tr=rep(c(0,1),each=sum(WW1)),
                  X=c(PF.fit[which(WW1),kk],PF.fit[which(WW1),kk]),M=1,distance.tolerance=0,ties=FALSE,Weight=2,
                  restrict=matrix(c(1:sum(WW1),(1:sum(WW1))+sum(WW1),rep(-1,sum(WW1))),nrow=sum(WW1),ncol=3,byrow=FALSE))
    Matchmat[which(WW1),thiscnames]<-matrix(c(which(WW1),which(WW1)[out111$index.control]),ncol=2,byrow=FALSE)

  }

  # row_num <- 0
  # # cname1<-c()
  # for(jj in 1:(trtnumber-1)){
  #   for(kk in (jj+1):trtnumber){
  #     row_num <- row_num+1
  #     # thistrt <- trtlevels[jj]
  #     # thattrt <- trtlevels[kk]
  #
  #     results_dfm$Trt1[row_num] <- trtlevels[jj]
  #     results_dfm$Trt2[row_num] <- trtlevels[kk]
  #     results_dfm$Param[row_num] <- nameContrast(trt1=results_dfm$Trt1[row_num], trt2=results_dfm$Trt2[row_num])
  #     results_dfm$Estimate[row_num] <- meanw[kk]-meanw[jj]
  #     results_dfm$Variance[row_num] <- (1/N)*(
  #       mean( (Yiw[,kk]-Yiw[,jj]-(results_dfm$Estimate[row_num]))^2 ) +
  #         mean( (Kiw^2+Kiw)*sigsqiw*
  #           (W==results_dfm$Trt1[row_num] | W==results_dfm$Trt2[row_num])
  #         )
  #     )
  #
  #     # varestimate[row_num]<-mean((Yiw[,kk]-Yiw[,jj]-(meanw[kk]-meanw[jj]))^2)+
  #     # mean((Kiw^2+Kiw)*sigsqiw*(W==thistrt | W==thattrt))
  #   }
  # }
  # if (row_num != taunumber) { stop("Error in for loop") }
  results_list <- estimateTau(
    trtlevels = trtlevels,
    meanw = meanw,
    trtnumber = trtnumber,
    taunumber = taunumber,
    N=N,
    #also get variance estimates
    Yiw=Yiw, Kiw=Kiw,sigsqiw=sigsqiw,W=W
  )

  tau_dfm <- results_list$tau_dfm
  # varestimate<-varestimate/N
  # names(tauestimate)<-cname1
  # names(varestimate)<-cname1
  # names(varestimateAI2012)<-cname1


  if (GPSM=="multinomiallogisticReg") {

    tau_dfm$VarianceAI2012 <- NA

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
            if(kkk==(kk+1)){C1.X1Y<-C1.X1Y*(1-PF.fit[,kk+1])}
            else if(kkk!=(kk+1))C1.X1Y<-C1.X1Y*(-PF.fit[,kk+1])
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
        result_row_num <- which(tau_dfm$Param  ==
                                  nameContrast(trt1=thistrt, trt2=thattrt))
        # cname1<-c(paste(paste(paste(paste(paste("EY(",thattrt,sep=""),")",sep=""),"-EY(",sep=""),thistrt,sep=""),")",sep=""))
        # varestimateAI2012[cname1]<-varestimate[cname1]-
        #   t(Cvec[jj,]+Cvec[kk,])%*%vcov_coeff%*%(Cvec[jj,]+Cvec[kk,])
        tau_dfm$VarianceAI2012[result_row_num] <-
          tau_dfm$Variance[result_row_num] -
            t(Cvec[jj,]+Cvec[kk,]) %*% vcov_coeff %*% (Cvec[jj,]+Cvec[kk,])
      }
    }
  }

  # tau_dfm <-
  #   as.data.frame(tau_dfm, stringsAsFactors = FALSE, row.names = NULL)
  # untidy_output <- list(tauestimate=tauestimate,
  #             varestimate=varestimate,
  #             varestimateAI2012=varestimateAI2012,
  #             analysis_idx=analysis_idx)

  # tidy_output <- tidyOutput(untidy_output=untidy_output)

  tidy_output <- list(
    results = tau_dfm,
    analysis_idx = analysis_idx,
    mu = results_list$mu_dfm,
    impute_mat = Yiw[prepared_data$sorted_to_orig,]
  )
  return(tidy_output)
}
