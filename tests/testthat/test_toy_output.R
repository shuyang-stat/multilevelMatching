library(multilevelMatching)
context("toy dataset results")


X <-c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9)
Y <-c(102,105,120,130,100,80,94,108,96)

W <- c(1, 1, 1, 3, 2, 3, 2, 1, 2)
existing_GPS_matrix <- cbind(
  c(0.5, 0.3, 0.5, 0.5, 0.5, 0.3, 0.3, 0.5, 0.3),
  c(1,1.6, 1, 1, 1, 1.6,1.6, 1,1.6)/6,
  c(2, 2.6, 2, 2, 2, 2.6, 2.6, 2, 2.6)/6
)

t1 <- multilevelMatchX(Y, W, X)
t2 <- multilevelGPSMatch(Y,W,X,Trimming=0,GPSM="multinomiallogisticReg")
t3 <- multilevelGPSMatch(Y,W,X,Trimming=1,GPSM="multinomiallogisticReg")
# t4 <- multilevelGPSMatch(Y=Y,W=W,X=existing_GPS_matrix,Trimming=0,GPSM="existing")
t_matX <- multilevelMatchX(Y, W, as.matrix(X))
# t_factorW <- multilevelMatchX(Y, as.factor(W), X)

my_tolerance <- 0.0001

test_that("match on one X returns same output", {
  expect_equal( (t1$results)$Estimate,
               c( -10.666667,   6.666667 , 17.333333),
               tolerance = my_tolerance)
  expect_equal( (t1$results)$Variance,
                c(  9.111111 ,615.580247, 613.925926),
                tolerance = my_tolerance)
})

test_that("match on GPS with one X and no trimming returns same output", {
  expect_equal( (t2$results)$Estimate,
                c(  -10.444444  , 6.666667 , 17.111111),
                tolerance = my_tolerance)
  expect_equal( (t2$results)$Variance,
                c( 8.545953, 616.913580 ,611.122085),
                tolerance = my_tolerance)
  expect_equal( (t2$results)$varestimateAI2012,
                c(   8.302024, 411.456234 ,434.247037),
                tolerance = my_tolerance)
})

test_that("match on GPS with one X with trimming returns same output", {
  expect_equal( (t3$results)$Estimate,
                c( -9.375 , 5.875, 15.250),
                tolerance = my_tolerance)
  expect_equal( (t3$results)$Variance,
                c(  7.794922 ,582.654297 ,576.304688),
                tolerance = my_tolerance)
  expect_equal( (t3$results)$varestimateAI2012,
                c(  5.072057 ,383.848575, 430.978089),
                tolerance = my_tolerance)
  # expect_equal( t3$analysis_idx,
  #               c(  5.072057 ,383.848575, 430.978089),
  #               tolerance = my_tolerance)
})
# test_that("match on GPS with existing GPS returns same output", {
#   expect_equal( (t4$results)$Estimate,
#                 c( -8.000000 , 1.777778 , 9.777778),
#                 tolerance = my_tolerance)
#   expect_equal( (t4$results)$Variance,
#                 c( 18.04938 ,573.25377, 552.78464),
#                 tolerance = my_tolerance)
# })
test_that("match on X with X a one-column matrix returns same output", {
  expect_equal( (t_matX$results)$Estimate,
                c( -10.666667 ,  6.666667 , 17.333333),
                tolerance = my_tolerance)
  expect_equal( (t_matX$results)$Variance,
                c( 9.111111 ,615.580247, 613.925926),
                tolerance = my_tolerance)
})
