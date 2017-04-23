library(multilevelMatching)
context("impute_mat is ordered as original data")


X <- matrix(c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9), ncol=1)
Y <- matrix(c(102,105,120,130,100,80,94,108,96), ncol=1)
W <- matrix(c(1, 1, 1, 3, 2, 3, 2, 1, 2), ncol=1)
GPSM <-   "multinomiallogisticReg"
Trimming <-  0

match_method <- "MatchOnGPS"

t2_in_imputemat <- multilevelGPSMatch(Y,W,X,Trimming=Trimming,GPSM=GPSM)

# prepared_data_dir <- file.path("..","prepared_toy_data.Rdata")
# load(prepared_data_dir)
# prepared_data_dir <- file.path("tests","prepared_toy_data.Rdata")
# ##GPS_imputes_mat <- t3$imputes_mat
# ### save(GPS_imputes_mat,file= tests_dir) #baseline_tests2 <- t2
# save(prepared_data, file= prepared_data_dir)

tests_dir <- file.path("..", "test_toy_output.Rdata")
# ##GPS_imputes_mat <- t3$imputes_mat
### save(GPS_imputes_mat,file= tests_dir)
load(tests_dir)


# existing_GPS_matrix <- cbind(
#   c(0.5, 0.3, 0.5, 0.5, 0.5, 0.3, 0.3, 0.5, 0.3),
#   c(1,1.6, 1, 1, 1, 1.6,1.6, 1,1.6)/6,
#   c(2, 2.6, 2, 2, 2, 2.6, 2.6, 2, 2.6)/6
# )

# t1 <- multilevelMatchX(Y, W, X)
# t3 <- multilevelGPSMatch(Y,W,X,Trimming=1,GPSM="multinomiallogisticReg")
# t4 <- multilevelGPSMatch(Y=Y,W=W,X=existing_GPS_matrix,Trimming=0,GPSM="existing")
# t_matX <- multilevelMatchX(Y, W, as.matrix(X))
# t_factorW <- multilevelMatchX(Y, as.factor(W), X)

# tests_dir <- file.path("..", "test_toy_output.Rdata")
# # tests_dir <- file.path("tests","test_toy_output.Rdata")
# # ##GPS_imputes_mat <- t3$imputes_mat
# ### save(GPS_imputes_mat,file= tests_dir) #baseline_tests2 <- t2
# # save(baseline_tests2, file= tests_dir)
# load(tests_dir)


my_tolerance <- 1e-3
test_that("match on GPS with one X and no trimming returns same output", {
  expect_equal( (t2_in_imputemat$results)$Estimate,
                (baseline_tests2$results)$Estimate,
                tolerance = my_tolerance)
  expect_equal( (t2_in_imputemat$results)$Variance,
                (baseline_tests2$results)$Variance,
                tolerance = my_tolerance)
  expect_equal( (t2_in_imputemat$results)$VarianceAI2012,
                (baseline_tests2$results)$VarianceAI2012,
                tolerance = my_tolerance)
  # expect_equal( (t2_in_imputemat$impute_mat),
  #               (baseline_tests2$impute_mat),
  #               tolerance = my_tolerance)
})


baseline_imputes <- matrix(NA,ncol=3,  nrow= length(W))
for (ii in 1:nrow(baseline_imputes)){
  baseline_imputes[ii,W[ii]] <- Y[ii]
}


new_imputes <- (t2_in_imputemat$impute_mat)
new_imputes[is.na(baseline_imputes)] <- NA

test_that("impute_mat behaves well", {
  expect_equal( new_imputes[1:3,],
                baseline_imputes[1:3,],
                tolerance = my_tolerance)
  expect_equal( new_imputes ,
                baseline_imputes ,
                tolerance = my_tolerance)
})


# load(prepared_data_dir)
# test_that("impute_mat behaves well, 2", {
# expect_equal( new_imputes[1:3,],
#               baseline_imputes[1:3,],
#               tolerance = my_tolerance)
# expect_equal( new_imputes ,
#               baseline_imputes ,
#               tolerance = my_tolerance)
# })
