

prepareData <- function(
  Y,W,X,match_method,
  Trimming = FALSE#Trimming_fit_args
){

  if (length(class(X))!=1 && class(X)!="matrix") {
    long_message <- "Dataframes, lists, and some other types are not supported at this time"
    warning("It is recommended that X is a matrix. X will be coerced to a matrix.")
    if (is.list(X) || is.data.frame(X)) {stop(long_message)}
    if (is.factor(X)) {
      warning("It is recommended that X is a matrix.\n Factor variables are not recommended at this time")
    }

    X <- try(as.matrix(X))
    if ("try-error" %in% class(X)){
      stop("Error with coercing X to matrix. Please specify X as a matrix")
    }

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
      analysisidx <- overlap.idx
    }

  } else {stop("match_method coding error")}

  # N <- length(Y) # number of observations ##moved into reorderByTreatment()

  ## order the treatment increasingly
  ordered_data_list <- reorderByTreatment(W=W,X=X,Y=Y)



  return(ordered_data_list)
}
