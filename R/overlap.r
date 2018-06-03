
#' Function to determine overlap from Crump
#'
#' @param PF.fit fitted propensity model
overlap<-function(PF.fit){
  obj<-function(alpha){
    gx<-apply(1/PF.fit,1,sum)
    id<-(gx<alpha)
    mean(gx*id)/(mean(id)^2)
  }
  gx<-apply(1/PF.fit,1,sum)
  #alpha1<-optim(median(gx), obj)$par
  suppressWarnings( alpha1<-stats::optim(stats::median(gx), obj)$par)
  alpha<-2*obj(alpha1)
  idx<-which(gx<alpha)
  return(list(alpha=alpha,idx=idx))
}
