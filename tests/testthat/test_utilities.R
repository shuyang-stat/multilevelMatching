
# context("new naming function")
context("Misc small improvements")

test_that(
  desc="nameCols() works",
  code = {

    trtnumber <- 3
    trtlevels <- c("foo", "var", "bar")
    cname<-c()
    for(kk in 1:trtnumber){
      thistrt<-trtlevels[kk]
      cname<-c(cname,c(paste(paste(paste("m",thistrt,sep=""),".",sep=""),1,sep=""),
                       paste(paste(paste("m",thistrt,sep=""),".",sep=""),2,sep="")))
    }

    expect_equal(
      cname,
      nameCols(trtlevels)
    )
  }
)




test_that(
  desc = "reorderByTreatment() can be made simpler + works with factors",
  {
    W <- 1:4
    temp1 <- sort(W, index.return=TRUE)
    temp2 <- list(x = W)
    temp2$ix <- 1:length(W)

    expect_equal(
      temp1,
      temp2
    )


    W <- (letters[1:5])
    temp1 <- sort(W, index.return=TRUE)
    temp2 <- list(x = W)
    temp2$ix <- 1:length(W)

    expect_equal(
      temp1,
      temp2
    )

    W <- as.factor(letters[1:5])
    temp1 <- sort(W, index.return=TRUE)
    expect_error(temp1$ix)
    ##use the below method!
    # temp2 <- list(x = W)
    # temp2$ix <- 1:length(W)
    #
    # expect_equal(
    #   temp1,
    #   temp2
    # )
  }
)
