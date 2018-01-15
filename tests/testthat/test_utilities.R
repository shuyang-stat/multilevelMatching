
context("new naming function")


trtnumber <- 3
trtlevels <- c("foo", "var", "bar")
cname<-c()
for(kk in 1:trtnumber){
  thistrt<-trtlevels[kk]
  cname<-c(cname,c(paste(paste(paste("m",thistrt,sep=""),".",sep=""),1,sep=""),
                   paste(paste(paste("m",thistrt,sep=""),".",sep=""),2,sep="")))
}



test_that(
  desc="nameCols() works",
  code = {
    expect_equal(
      cname,
      nameCols(trtlevels)
    )
  }
)
