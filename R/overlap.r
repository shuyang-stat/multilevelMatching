
# #' Function to determine overlap from Crump et al. (2009)'s method.
# #'
# #' @param PF.fit fitted propensity model
# #'
# #' @references Crump, R. K., Hotz, V. J., Imbens, G. W., & Mitnik, O. A. (2009).
# #'   Dealing with limited overlap in estimation of average treatment effects.
# #'   Biometrika, 96(1), 187-199. \url{https://doi.org/10.1093/biomet/asn055}
# #'
overlap <- function(PF.fit){

  obj <- function(alpha){
    gx <- apply(1/PF.fit,1,sum)
    id <- (gx<alpha)
    mean(gx*id)/(mean(id)^2)
  }
  gx <- apply(1/PF.fit,1,sum)
  #alpha1<-optim(median(gx), obj)$par
  suppressWarnings( alpha1 <- stats::optim(stats::median(gx), obj)$par)
  alpha <- 2*obj(alpha1)
  idx <- which(gx<alpha)

  return(list(alpha=alpha,idx=idx))
}
