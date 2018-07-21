#' Matching on X with multilevel treatments
#'
#' @param Y A continuous response vector (1 x n)
#' @param W A treatment vector (1 x n) with numerical values indicating
#'   treatment groups
#' @param X A covariate matrix (p x n) with no intercept
#'
#' @return A list with 2 elements: \code{tauestimate}, \code{varestimate}, where
#'   \code{tauestimate} is a vector of estimates for pairwise treatment effects,
#'   and \code{varestimate} is a vector of variance estimates for
#'   \code{tauestimate}, using Abadie & Imbens (2006)'s method.
#'
#' @references Yang, S., Imbens G. W., Cui, Z., Faries, D. E., & Kadziola, Z.
#'   (2016) Propensity Score Matching and Subclassification in Observational
#'   Studies with Multi-Level Treatments. Biometrics, 72, 1055-1065.
#'   \url{https://doi.org/10.1111/biom.12505}
#'
#'   Abadie, A., & Imbens, G. W. (2006). Large sample properties of matching
#'   estimators for average treatment effects. econometrica, 74(1), 235-267.
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
#' @seealso \code{\link{multilevelGPSMatch}};
#'   \code{\link{multilevelGPSStratification}}
#'
#' @examples
#'   X<-c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9)
#'   Y<-c(102,105,120,130,100,80,94,108,96)
#'   W<-c(1,1,1,3,2,3,2,1,2)
#'   multilevelMatchX(Y,W,X)
#'
#' @export
multilevelMatchX <- function(Y,W,X){

  ## order the treatment increasingly
  if(1-is.unsorted(W)){
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

  # number of treatment levels
  trtnumber <- length(unique(W))
  # all treatment levels
  trtlevels <- unique(W)
  # number of observations by treatment level
  pertrtlevelnumber <- table(W)
  # number of pairwise treatment effects
  taunumber <- (trtnumber*(trtnumber+1)/2)-trtnumber

  # prepared_data <- prepareData_legacy(
  #   Y=Y, W=W, X=X,
  #   match_method = match_method#,
  #   # Trimming = FALSE
  #   #Trimming_fit_args
  # )
  # W <- prepared_data$W
  # X <- prepared_data$X
  # Y <- prepared_data$Y
  # N <- prepared_data$N
  # trtnumber <- prepared_data$trtnumber
  # trtlevels <- prepared_data$trtlevels
  # pertrtlevelnumber <- prepared_data$pertrtlevelnumber
  # taunumber <- prepared_data$taunumber
  # analysis_idx <- prepared_data$analysis_idx

  tauestimate <- varestimate <- rep(NA,taunumber)
  meanw <- rep(NA,trtnumber)
  ## Yiw is the full imputed data set
  Yiw <- matrix(NA,N,trtnumber)
  ## Kiw is vector of number of times unit i used as a match
  Kiw <- sigsqiw <- matrix(NA,N,1)

  Matchmat <- matrix(NA,N,trtnumber*2)
  cname <- c()
  for(kk in 1:trtnumber){
    thistrt <- trtlevels[kk]
    cname <- c(cname,c(paste(paste(paste("m",thistrt,sep=""),".",sep=""),1,sep=""),
                     paste(paste(paste("m",thistrt,sep=""),".",sep=""),2,sep="")))
  }
  colnames(Matchmat) <- cname


  for(kk in 1:trtnumber){

    thistrt <- trtlevels[kk]
    if(kk==1){fromto <- 1:pertrtlevelnumber[1]}
    if(kk>1){fromto <- (1:pertrtlevelnumber[kk])+sum(pertrtlevelnumber[1:(kk-1)])}
    W1 <- W!=thistrt
    out1 <- Matching::Match(Y=Y,Tr=W1,X=X,distance.tolerance=0,ties=FALSE,Weight=2)
    mdata1 <- out1$mdata
    meanw[kk] <- stats::weighted.mean(c(Y[which(W==thistrt)],mdata1$Y[which(mdata1$Tr==0)]),c(rep(1,length(which(W==thistrt))),out1$weights))
    Kiw[fromto,1] <- table(factor(out1$index.control,levels=fromto))
    Yiw[which(W==thistrt),kk] <- Y[which(W==thistrt)]
    Yiw[which(W!=thistrt),kk] <- mdata1$Y[which(mdata1$Tr==0)]
    WW1 <- W==thistrt
    X <- as.matrix(X)
    out11 <- Matching::Match(Y=rep(Y[which(WW1)],times=2),Tr=rep(c(1,0),each=sum(WW1)),
                 X=rbind(as.matrix(X[which(WW1),]),as.matrix(X[which(WW1),])),M=1,distance.tolerance=0,ties=FALSE,Weight=2,
                 restrict=matrix(c(1:sum(WW1),(1:sum(WW1))+sum(WW1),rep(-1,sum(WW1))),nrow=sum(WW1),ncol=3,byrow=FALSE))

    mdata11 <- out11$mdata
    temp11 <- (mdata11$Y[which(mdata11$Tr==1)]-mdata11$Y[which(mdata11$Tr==0)])^2/2
    sigsqiw[which(W==thistrt),1] <- temp11

    thiscnames <- c(paste(paste(paste("m",thistrt,sep=""),".",sep=""),1,sep=""),
                  paste(paste(paste("m",thistrt,sep=""),".",sep=""),2,sep=""))

    # find two outsiders closest
    findmatch1 <- Matching::Match(Y=Y,Tr=W1,X=X,distance.tolerance=0,ties=FALSE,Weight=2,M=2)
    Matchmat[unique(findmatch1$index.treated),thiscnames] <- matrix(findmatch1$index.control,ncol=2,byrow=TRUE)
    # find one insider closest
    out111 <- Matching::Match(Y=rep(Y[which(WW1)],times=2),Tr=rep(c(0,1),each=sum(WW1)),
                  X=rbind(as.matrix(X[which(WW1),]),as.matrix(X[which(WW1),])),M=1,distance.tolerance=0,ties=FALSE,Weight=2,
                  restrict=matrix(c(1:sum(WW1),(1:sum(WW1))+sum(WW1),rep(-1,sum(WW1))),nrow=sum(WW1),ncol=3,byrow=FALSE))
    Matchmat[which(WW1),thiscnames] <- matrix(c(which(WW1),which(WW1)[out111$index.control]),ncol=2,byrow=FALSE)

  }

  cnt <- 0
  cname1 <- c()
  for(jj in 1:(trtnumber-1)){
    for(kk in (jj+1):trtnumber){
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

  out_list <- list(
    tauestimate = tauestimate,
    varestimate = varestimate
  )
  return(out_list)


  # estimate_args <- list(
  #   trtlevels = trtlevels,
  #   meanw = meanw,
  #   trtnumber = trtnumber,
  #   taunumber = taunumber,
  #   N=N,
  #   #also get variance estimates
  #   Yiw=Yiw, Kiw=Kiw,sigsqiw=sigsqiw,W=W
  # )
  # results_list <- do.call(estimateTau_legacy,estimate_args)
  #
  #
  # tau_dfm <- results_list$tau_dfm
  #
  #
  # # untidy_output <- list(tauestimate=tauestimate,varestimate=varestimate)
  # # tidy_output <- tidyOutput(untidy_output=untidy_output)
  #
  #
  # tidy_output <- list(
  #   results = tau_dfm,
  #   analysis_idx = analysis_idx,
  #   mu = results_list$mu_dfm,
  #   impute_mat = Yiw[prepared_data$sorted_to_orig,],
  #   estimate_args = estimate_args
  # )
  # return(tidy_output)
}

