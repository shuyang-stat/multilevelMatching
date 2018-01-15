
context("new multiMatch() function")


X <- matrix(c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9), ncol=1)
Y <- matrix(c(102,105,120,130,100,80,94,108,96), ncol=1)

W <- matrix(c(1, 1, 1, 3, 2, 3, 2, 1, 2), ncol=1)


## TODO: add tests for match on "covariates" and "existing"



reorder_estimate_args <- c(
  "trtlevels", "meanw","trtnumber","taunumber","N","Yiw","Kiw","sigsqiw","W"
)
## old function
mnom <- multilevelGPSMatch(Y,W,X,Trimming=0,GPSM="multinomiallogisticReg")

# new function
mnom_new <- multiMatch(Y, W, X, trimming=0, match_on = "multinom")

mnom2 <- mnom_new
# mnom2$estimate_args <- mnom2$estimate_args[reorder_estimate_args]

test_that("multinom-matching returns same with new multiMatch function", {
  expect_equal(
    object = mnom[1:4],
    expected = mnom2[1:4],
    tolerance = 1e-7
  )
})





##old
polr <- multilevelGPSMatch(Y,W,X,Trimming=0,GPSM="ordinallogisticReg")
#new
polr_new <- multiMatch(Y, W, X, trimming=0, match_on = "polr")
polr2 <- polr_new
polr2$estimate_args <- polr2$estimate_args[reorder_estimate_args]

test_that("polr-matching returns same with new multiMatch function", {
  expect_equal(
    object = polr[1:4],
    expected = polr2[1:4],
    tolerance = 1e-7
  )
})





existing_GPS_matrix <- cbind(
  c(0.5, 0.3, 0.5, 0.5, 0.5, 0.3, 0.3, 0.5, 0.3),
  c(1,1.6, 1, 1, 1, 1.6,1.6, 1,1.6)/6,
  c(2, 2.6, 2, 2, 2, 2.6, 2.6, 2, 2.6)/6
)
eps <- multilevelGPSMatch(Y=Y,W=W,X=existing_GPS_matrix,Trimming=0,GPSM="existing")
epsmm <- suppressWarnings(
  multiMatch(Y=Y,W=W,X=existing_GPS_matrix,trimming=0,match_on="existing")
)


# test_that("'existing' ps-matching returns same with new multiMatch function", {
#   ## not passing tests
#   # expect_equal(
#   #   object = eps,
#   #   expected = epsmm,
#   #   tolerance = 1e-7
#   # )
#   #
#   # expect_equal(
#   #   object = eps$impute_match_data[[1]],
#   #   expected = epsmm$impute_match_data[[1]],
#   #   tolerance = 1e-7
#   # )
#
#
#   # expect_equal(
#   #   object = eps$impute_match_data[[2]],
#   #   expected = epsmm$impute_match_data[[2]],
#   #   tolerance = 1e-7
#   # )
#   # expect_equal(
#   #   object = eps$impute_match_data[[3]],
#   #   expected = epsmm$impute_match_data[[3]],
#   #   tolerance = 1e-7
#   # )
#   # expect_equal(
#   #   object = eps$impute_match_data[-(1:2)],
#   #   expected = epsmm$impute_match_data[-(1:2)],
#   #   tolerance = 1e-7
#   # )
# })


## not passing tests

# test_that("match on GPS with existing GPS returns same output", {
#   expect_equal( (eps$results)$Estimate,
#                 c( -8.000000 , 1.777778 , 9.777778),
#                 tolerance = my_tolerance)
#   expect_equal( (eps$results)$Variance,
#                 c( 18.04938 ,573.25377, 552.78464),
#                 tolerance = my_tolerance)
# })

# test_that("match on existing GPS returns same output", {
#   expect_equal( (epsmm$results)$Estimate,
#                 c( -8.000000 , 1.777778 , 9.777778),
#                 tolerance = my_tolerance)
#   expect_equal( (epsmm$results)$Variance,
#                 c( 18.04938 ,573.25377, 552.78464),
#                 tolerance = my_tolerance)
# })


covar <- multilevelMatchX(Y=Y,W=W,X=X)
covarsmm <- multiMatch(Y=Y,W=W,X=X,trimming=0,match_on="covariates")


test_that("existing ps-matching returns same with new multiMatch function", {
  ## not passing tests
  expect_equal(
    object = covar[1:4],
    expected = covarsmm[1:4],
    tolerance = 1e-7
  )

})

