#' Matching on X with multilevel treatments
#'
#' @param Y a continuous response vector (1 x n)
#' @param W a treatment vector (1 x n) with numerical values indicating treatment groups
#' @param X a covariate matrix (p x n) with no intercept
#'
#' @return according to \code{\link{estimateTau}}: a dataframe with two columns,
#' tauestimate, varestimate, where
#' tauestimate is a vector of estimates for pairwise treatment effects, and
#' varestimate is a vector of variance estimates for tauestimate, using
#' Abadie&Imbens(2006)'s method.
#'
#' @seealso \code{\link{multilevelGPSMatch}}; \code{\link{multilevelGPSStratification}}
#'
#' @examples
#'   X<-c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9)
#'   Y<-c(102,105,120,130,100,80,94,108,96)
#'   W<-c(1,1,1,3,2,3,2,1,2)
#'   multilevelMatchX(Y,W,X)
#'
#' @import Matching boot nnet MASS
#'
#' @export
multilevelMatchX <- function(Y,W,X){

  ## some checks
  match_method <- "MatchOnX"

  prepared_data <- prepareData(
    Y=Y, W=W, X=X,
    match_method = match_method#,
    # Trimming = FALSE
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

  # tauestimate<-varestimate<-rep(NA,taunumber)
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

  # cnt<-0
  # cname1<-c()
  # for(jj in 1:(trtnumber-1)){
  #   for(kk in (jj+1):trtnumber){
  #     cnt<-cnt+1
  #     thistrt<-trtlevels[jj]
  #     thattrt<-trtlevels[kk]
  #     cname1<-c(cname1,paste(paste(paste(paste(paste("EY(",thattrt,sep=""),")",sep=""),"-EY(",sep=""),thistrt,sep=""),")",sep=""))
  #     tauestimate[cnt]<-meanw[kk]-meanw[jj]
  #     varestimate[cnt]<-mean((Yiw[,kk]-Yiw[,jj]-(meanw[kk]-meanw[jj]))^2)+mean((Kiw^2+Kiw)*sigsqiw*(W==thistrt | W==thattrt))
  #   }
  # }
  # varestimate<-varestimate/N
  # names(tauestimate)<-cname1
  # names(varestimate)<-cname1

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


  # untidy_output <- list(tauestimate=tauestimate,varestimate=varestimate)
  # tidy_output <- tidyOutput(untidy_output=untidy_output)


  tidy_output <- list(
    results = tau_dfm,
    analysis_idx = analysis_idx,
    mu = results_list$mu_dfm,
    impute_mat = Yiw[prepared_data$sorted_to_orig,]
  )
  return(tidy_output)
}

