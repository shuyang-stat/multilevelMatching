#' Stratification on GPS with multilevel treatments
#'
#' @param Y a continuous response vector (1 x n)
#' @param W a treatment vector (1 x n) with numerical values indicating treatment groups
#' @param X a covariate matrix (p x n) with no intercept
#' @param GPSM an indicator of the methods used for estimating GPS, options include "multinomiallogisticReg", "ordinallogisticReg", and "existing"
#' @param NS (only required in the function multilevelGPSStratification) the number of strata
#' @param linearp (only required in the function multilevelGPSStratification) an indicator of subclassification on GPS (=0) or linear predictor of GPS (=1)
#' @param nboot (only required in the function multilevelGPSStratification) the number of boot replicates for variance estimation
#' @param model_options A list of the options to pass to propensity model.
#'   Currently under development. Can only pass reference level to multinomial
#'   logisitc regression.
#'
#' @return according to \code{\link{estimateTau}}: a dataframe with two columns,
#' tauestimate, varestimate, where
#' tauestimate is a vector of estimates for pairwise treatment effects, and
#' varestimate is a vector of variance estimates, using bootstrapping method.
#'
#' @seealso \code{\link{multilevelGPSMatch}}; \code{\link{multilevelMatchX}}
#'
#'
#' @import  boot nnet MASS
#'
#' @export
multilevelGPSStratification <- function(
  Y,W,X,NS,GPSM="multinomiallogisticReg",linearp=0,nboot,
  model_options = list(reference_level = "1") ){
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

  prepared_data <- prepareData_legacy(
    Y=Y, W=W, X=X,
    match_method = match_method#,
    # Trimming = Trimming
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
    if (linearp==1) {
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
    ## bug checks migrated to prepareData()
    PF.fit <- X
  }


  # results_dfm <- list(
  #   Param = rep(NA, taunumber),
  #   Trt1 = rep(NA, taunumber),
  #   Trt2 = rep(NA, taunumber),
  #   Estimate = rep(NA, taunumber),
  #   Variance = rep(NA, taunumber)#,
  #   # VarianceAI2012 = rep(NA,taunumber)
  # )

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

  # tauestimate<-rep(NA,taunumber)

  # cnt<-0
  # cname1<-c()

  ## Calculate estimates for tau
  results_list <- estimateTau_legacy(
    trtlevels = trtlevels,
    meanw = meanw,
    trtnumber = trtnumber,
    taunumber = taunumber,
    N=N
    # do NOT get variance estimates in stratification
    # Yiw=Yiw, Kiw=Kiw,sigsqiw=sigsqiw,W=W
  )
  tau_dfm <- results_list$tau_dfm
  # row_num <- 0
  # for(jj in 1:(trtnumber-1)){
  #   for(kk in (jj+1):trtnumber){
  #     # cnt<-cnt+1
  #     # thistrt<-trtlevels[jj]
  #     # thattrt<-trtlevels[kk]
  #     # cname1<-c(cname1,paste(paste(paste(paste(paste("EY(",thattrt,sep=""),")",sep=""),"-EY(",sep=""),thistrt,sep=""),")",sep=""))
  #     # tauestimate[cnt]<-meanw[kk]-meanw[jj]
  #     row_num <- row_num+1
  #     results_dfm$Trt1[row_num] <- trtlevels[jj]
  #     results_dfm$Trt2[row_num] <- trtlevels[kk]
  #     results_dfm$Param[row_num] <- nameContrast(trt1=results_dfm$Trt1[row_num], trt2=results_dfm$Trt2[row_num])
  #     results_dfm$Estimate[row_num] <- meanw[kk]-meanw[jj]
  #   }
  # }
  # # names(tauestimate)<-cname1
  # if (row_num != taunumber) { stop("Error in for loop") }


  ## Bootstrap the variance
  data<-cbind(W,Y,X)
  results <- boot::boot(data=data, statistic=estforboot,R=nboot,
                  GPSM=GPSM,linearp=linearp,trtnumber=trtnumber,
                  trtlevels=trtlevels,taunumber=taunumber,NS=NS)
  bootvar <- apply(results$t,2,stats::var,na.rm = TRUE)
  # names(bootvar)<-cname1

  ## Tidy the output
  tau_dfm$Variance <- bootvar
  # untidy_output <- list(tauestimate=tauestimate,varestimate=bootvar)
  # tidy_output <- tidyOutput(untidy_output=untidy_output)
  tau_dfm <-
    as.data.frame(tau_dfm, stringsAsFactors = FALSE, row.names = NULL)


  tidy_output <- list(
    results = tau_dfm,
    analysis_idx = analysis_idx,
    mu = results_list$mu_dfm
  )
  return(tidy_output)
}
