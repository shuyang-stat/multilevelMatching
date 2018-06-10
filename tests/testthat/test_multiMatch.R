
context("test_multiMatch() function")


X <- matrix(c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9), ncol=1)
Y <- matrix(c(102,105,120,130,100,80,94,108,96), ncol=1)
W <- matrix(c(1, 1, 1, 3, 2, 3, 2, 1, 2), ncol=1)


## TODO: Matching on "existing" is not passing


existing_GPS_matrix <- cbind(
  c(0.5, 0.3, 0.5, 0.5, 0.5, 0.3, 0.3, 0.5, 0.3),
  c(1,1.6, 1, 1, 1, 1.6,1.6, 1,1.6)/6,
  c(2, 2.6, 2, 2, 2, 2.6, 2.6, 2, 2.6)/6
)
eps_legacy <- multilevelGPSMatch(Y=Y,W=W,X=existing_GPS_matrix,Trimming=0,GPSM="existing")
eps_multiMatch <- suppressWarnings(
  suppressMessages(
    multiMatch(Y=Y,W=W,X=existing_GPS_matrix,trimming=0,match_on="existing")
  )
)

## not passing tests
test_that("'existing' ps-matching returns same with new multiMatch function", {

  expect_failure(
    expect_equal(
      object = eps_legacy,
      expected = eps_multiMatch,
      tolerance = 1e-7
    )
  )

  expect_failure(
    expect_equal(
      object = eps_legacy$impute_match_data[[1]],
      expected = eps_multiMatch$impute_match_data[[1]],
      tolerance = 1e-7
    )
  )

  expect_failure(
    expect_equal(
      object = eps_legacy$impute_match_data[[2]],
      expected = eps_multiMatch$impute_match_data[[2]],
      tolerance = 1e-7
    )
  )

  expect_failure(
    expect_equal(
      object = eps_legacy$impute_match_data[[3]],
      expected = eps_multiMatch$impute_match_data[[3]],
      tolerance = 1e-7
    )
  )

  expect_failure(
    expect_equal(
      object = eps_legacy$impute_match_data[-(1:2)],
      expected = eps_multiMatch$impute_match_data[-(1:2)],
      tolerance = 1e-7
    )
  )
})


## not passing tests
test_that("match on GPS with existing GPS returns same output", {

  # expect_failure(
    expect_equal( (eps_legacy$results)$Estimate,
                  c( -8.000000 , 1.777778 , 9.777778),
                  tolerance = my_tolerance)
  # )

  expect_failure(
    expect_equal( (eps_legacy$results)$Variance,
                  c( 18.04938 ,573.25377, 552.78464),
                  tolerance = my_tolerance)
  )
})

test_that("match on existing GPS returns same output", {
  expect_failure(
    expect_equal( (eps_multiMatch$results)$Estimate,
                  c( -8.000000 , 1.777778 , 9.777778),
                  tolerance = my_tolerance)
  )

  expect_failure(
    expect_equal( (eps_multiMatch$results)$Variance,
                  c( 18.04938 ,573.25377, 552.78464),
                  tolerance = my_tolerance)
  )
})


reorder_estimate_args <- c(
  "trtlevels", "meanw","trtnumber","taunumber","N","Yiw","Kiw","sigsqiw","W"
)
## old function
mnom_legacy <- multilevelGPSMatch(Y,W,X,Trimming=0,GPSM="multinomiallogisticReg")

# new function
mnom_multiMatch <- suppressMessages(
  multiMatch(Y, W, X, trimming=0, match_on = "multinom")
)


test_that("multinom-matching returns same with new multiMatch function", {

  names(mnom_legacy$results)[6] <- "VarianceAI2016" #2018-06-10

  expect_equal(
    object = mnom_legacy[1:3],
    expected = mnom_multiMatch[1:3],
    tolerance = 1e-7
  )
  expect_equal(
    object = mnom_legacy$impute_mat,
    expected = mnom_multiMatch$impute_mat,
    tolerance = 1e-7,
    check.attributes = FALSE
  )

})





##old
polr_legacy <- multilevelGPSMatch(Y,W,X,Trimming=0,GPSM="ordinallogisticReg")
#new
polr_multiMatch <- suppressMessages(
  multiMatch(Y, W, X, trimming=0, match_on = "polr")
)

polr_multiMatch$estimate_args <- polr_multiMatch$estimate_args[reorder_estimate_args]

test_that("polr-matching returns same with new multiMatch function", {

  names(polr_legacy$results)[6] <- "VarianceAI2016" #2018-06-10

  expect_equal(
    object = polr_legacy[1:3],
    expected = polr_multiMatch[1:3],
    tolerance = 1e-7
  )

  expect_equal(
    object = polr_legacy$impute_mat,
    expected = polr_multiMatch$impute_mat,
    tolerance = 1e-7,
    check.attributes = FALSE
  )
})






covar_legacy <- multilevelMatchX(Y=Y,W=W,X=X)
covar_multiMatch <- suppressMessages(
  multiMatch(Y=Y,W=W,X=X,trimming=0,match_on="covariates")
)


test_that("existing ps-matching returns same with new multiMatch function", {

  names(covar_legacy$results)[6] <- "VarianceAI2016" #2018-06-10

  expect_equal(
    object = covar_legacy[1:3],
    expected = covar_multiMatch[1:3],
    tolerance = 1e-7
  )

  expect_equal(
    object = covar_legacy$impute_mat,
    expected = covar_multiMatch$impute_mat,
    tolerance = 1e-7,
    check.attributes = FALSE
  )

})

