

prepareData <- function(
  Y,W,X,match_method,
  Trimming = FALSE#Trimming_fit_args
){

  if (class(X)!="matrix"){
    stop("X must be of class 'matrix'.\n
         Dataframes, lists, and other types are not supported at this time")
  }
  ## some checks
  argChecks(Y=Y, W=W, X=X, N=NULL)


  if (Trimming==0) { analysisidx <- 1:length(Y) }

  if (match_method == "MatchOnX") {
    #nothing
  } else if (match_method == "StratifyOnGPS") {
    #nothing
  } else if (match_method == "MatchOnGPS"){

    if(Trimming==1){
      #PF modeling
      W.ref <- relevel(as.factor(W),ref=1)
      temp<-capture.output(PF.out <- multinom(W.ref~X))  ##make this into a subfunction?
      PF.fit <- fitted(PF.out)
      ## identify sufficient overlap
      overlap.idx<-overlap(PF.fit)$idx
      W <- W[overlap.idx]
      X <- X[overlap.idx,]
      Y <- Y[overlap.idx]
      analysisidx<-overlap.idx
    }

  } else {stop("match_method coding error")}

  # N <- length(Y) # number of observations ##moved into reorderByTreatment()

  ## order the treatment increasingly
  ordered_data_list <- reorderByTreatment(W=W,X=X,Y=Y)



  return(ordered_data_list)
}
