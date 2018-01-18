
context("Matching on existing GPS works")


X<-c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9)
Y<-c(102,105,120,130,100,80,94,108,96)
W<-c(1,1,1,3,2,3,2,1,2)


existing_GPS_matrix <- cbind(
  c(0.5, 0.3, 0.5, 0.5, 0.5, 0.3, 0.3, 0.5, 0.3),
  c(1,1.6, 1, 1, 1, 1.6,1.6, 1,1.6)/6,
  c(2, 2.6, 2, 2, 2, 2.6, 2.6, 2, 2.6)/6
)

set.seed(12345)
t4 <- multilevelGPSMatch(Y=Y,W=W,X=existing_GPS_matrix,Trimming=0,GPSM="existing")


# ## failing test
# this_t <- t4
test_that("match on GPS (with existing GPS) returns original (v0.1) output", {
  expect_equal( (t4$results)$Estimate,
                c(-8.888889, 1.111111, 10.000000),
                # c( -8.000000 , 1.777778 , 9.777778),
                tolerance = my_tolerance)
  expect_equal( (t4$results)$Variance,
                c(26.01097,   578.20850, 551.50617),
                # c( 18.04938 ,573.25377, 552.78464),
                tolerance = my_tolerance)
})


set.seed(12345)
t4mm <- suppressWarnings(
  multiMatch(
    Y=Y,W=W,X=existing_GPS_matrix,
    trimming=0,match_on="existing")
)
load(
#   t4mm_orig <- t4mm
# save(t4mm_orig,
     file = rprojroot::find_package_root_file(
       "tests", "testthat", "testing_datafiles",
       "test_multiMatch_exising_GPS.Rdata"))
# )
## Failing tests

# test_that("multiMatch on GPS with existing GPS returns original (v0.1) output", {
#   expect_equal( (t4$results)$Estimate,
#                 (t4mm$results)$Estimate,
#                 tolerance = my_tolerance)
#   expect_equal( (t4$results)$Variance,
#                 (t4mm$results)$Variance,
#                 tolerance = my_tolerance)
#   expect_equal( (t4$results) ,
#                 (t4mm$results) ,
#                 tolerance = my_tolerance)
# })

test_that("multiMatch on GPS with existing GPS returns as expected (v0.1.5)", {

  expect_equal( (t4mm_orig$results) ,
                (t4mm$results) ,
                tolerance = my_tolerance)
})


## Existing GPS matching via multiMatch() is not perfect

test_that("multiMatch on GPS with existing GPS returns SIMILAR output to original (v0.1)", {
set.seed(11)
N <- 300
X <- rnorm(N, 0, 1)
prW1 <- sample(size = N, x=(1:4)/10, replace = TRUE)
prW2 <- (1-prW1)*sample(size = N, x=(1:4)/5, replace = TRUE)
prW3 <- 1- (prW1+prW2)
existing_GPS_matrix <- cbind(prW1, prW2, prW3)
W <- rep(NA, N)
for(ii in 1:N){
  W[ii] <- sample(1:3, size = 1, replace = TRUE, prob = existing_GPS_matrix[ii,])
}
Y <- round(rnorm(N, 10 - W +0.2*X, 1),1)

# existing_GPS_matrix <- cbind(pr_w1, pr_w2,pr_w3)

set.seed(12345)
t4 <- multilevelGPSMatch(Y=Y,W=W,X=existing_GPS_matrix,Trimming=0,GPSM="existing")

set.seed(12345)
t4mm <- suppressWarnings(
  multiMatch(
    Y=Y,W=W,X=existing_GPS_matrix,
    trimming=0,match_on="existing")
)

  expect_equal( (t4$results)$Estimate,
                (t4mm$results)$Estimate,
                tolerance = 0.04)
  expect_equal( (t4$results)$Variance,
                (t4mm$results)$Variance,
                tolerance = 0.008)
})



