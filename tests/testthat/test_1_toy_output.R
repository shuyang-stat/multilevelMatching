
context("toy dataset results")


X <- matrix(c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9), ncol=1)
Y <- matrix(c(102,105,120,130,100,80,94,108,96), ncol=1)
rownames(Y) <- letters[4+(1:9)]
W <- matrix(c(1, 1, 1, 3, 2, 3, 2, 1, 2), ncol=1)
existing_GPS_matrix <- cbind(
  c(0.5, 0.3, 0.5, 0.5, 0.5, 0.3, 0.3, 0.5, 0.3),
  c(1,1.6, 1, 1, 1, 1.6,1.6, 1,1.6)/6,
  c(2, 2.6, 2, 2, 2, 2.6, 2.6, 2, 2.6)/6
)

my_tolerance <- 0.0001
Param_names <- c( "EY(2)-EY(1)", "EY(3)-EY(1)" ,"EY(3)-EY(2)")
Trt1s <- c(1,1,2)
Trt2s <- c(2,3,3)



# t4 <- multilevelGPSMatch(Y=Y,W=W,X=existing_GPS_matrix,Trimming=0,GPSM="existing")
## Tests for matching on existing GPS are
##    moved to test_existing_GPS_matching.R
##
## In summary, multiMatch() can sometimes produce slightly different results
##    than multilevelGPSMatching(), like when there are ties.


# t_factorW <- multilevelMatchX(Y, as.factor(W), X)
## Tests for when W is a factor not yet implemented; issue BarkleyBG/#1


test_that("multilevelMatchX on one X returns same output as v0.1", {
  t1 <- multilevelMatchX(Y, W, X)
  this_t <- t1

  expect_equal(
    object = (this_t$results)$Estimate,
    expected = c( -10.666667,   6.666667 , 17.333333),
    tolerance = my_tolerance)
  expect_equal(
    object =(this_t$results)$Variance,
    expected = c(  9.111111 ,615.580247, 613.925926),
    tolerance = my_tolerance)
  expect_identical(
    object = (this_t$results)$Param, expected = Param_names)
  expect_identical( (this_t$results)$Trt1, Trt1s)
  expect_identical( (this_t$results)$Trt2, Trt2s)
})

## More tests between multiMatch and multilevelMatchX in other files
test_that("multiMatch returns same as multilevelMatchX on one X", {

  t1 <- multilevelMatchX(Y, W, X)
  t1_v2 <- multiMatch(Y, W, X, match_on = "covariates")

  expect_equal(
    object = t1[1:3],
    expected = t1_v2[1:3],
    tolerance = 1e-7
  )
  expect_equal(
    object = t1$impute_mat,
    expected = t1_v2$impute_mat,
    tolerance = 1e-7,
    check.attributes = FALSE
  )
})


# reorder_estimate_args <- c(
#   "trtlevels", "meanw","trtnumber","taunumber","N","Yiw","Kiw","sigsqiw","W"
# )




test_that(
  "multilevelGPSMatch with one X, no tri, returns same output as v0.1",
  {
    t2 <- multilevelGPSMatch(Y,W,X,Trimming=0,GPSM="multinomiallogisticReg")
    this_t <- t2

    tests_data <- quickLookup("test_toy_output.Rdata")
    load(tests_data)



    expect_equal( (this_t$results)$Estimate,
                  c(  -10.444444  , 6.666667 , 17.111111),
                  tolerance = my_tolerance)
    expect_equal( (this_t$results)$Variance,
                  c( 8.545953, 616.913580 ,611.122085),
                  tolerance = my_tolerance)
    expect_equal( (this_t$results)$VarianceAI2012,
                  c(   8.302024, 411.456234 ,434.247037),
                  tolerance = my_tolerance)
    expect_identical( (this_t$results)$Param, Param_names)
    expect_identical( (this_t$results)$Trt1, Trt1s)
    expect_identical( (this_t$results)$Trt2, Trt2s)
    expect_equal(
      object = (this_t$impute_mat),
      expected = baseline_tests2$impute_mat,
      tolerance = my_tolerance
    )
  }
)

test_that(
  "multilevelGPSMatch one X and trimming returns same output as v0.1",
  {
    t3 <- multilevelGPSMatch(Y,W,X,Trimming=1,GPSM="multinomiallogisticReg")
    this_t <- t3


    expect_equal( (this_t$results)$Estimate,
                  c( -9.375 , 5.875, 15.250),
                  tolerance = my_tolerance)
    expect_equal(
      object = (this_t$results)$Variance,
      expected =  c(  7.794922 ,582.654297 ,576.304688),
      tolerance = my_tolerance)
    expect_equal( (this_t$results)$VarianceAI2012,
                  c(  5.072057 ,383.848575, 430.978089),
                  tolerance = my_tolerance)
    expect_identical( (this_t$results)$Param, Param_names)
    expect_identical( (this_t$results)$Trt1, Trt1s)
    expect_identical( (this_t$results)$Trt2, Trt2s)

    # expect_equal( this_t$analysis_idx,
    #               c(  5.072057 ,383.848575, 430.978089),
    #               tolerance = my_tolerance)
  }
)


test_that(
  "multilevelMatchX with one-column matrix X returns same output as v0.1",
  {
    t_matX <- multilevelMatchX(Y, W, as.matrix(X))
    this_t <- t_matX

    expect_equal(
      (this_t$results)$Estimate,
      c( -10.666667 ,  6.666667 , 17.333333),
      tolerance = my_tolerance
    )
    expect_equal(
      (this_t$results)$Variance,
      c( 9.111111 ,615.580247, 613.925926),
      tolerance = my_tolerance
    )
    expect_identical( (this_t$results)$Param, Param_names)
    expect_identical( (this_t$results)$Trt1, Trt1s)
    expect_identical( (this_t$results)$Trt2, Trt2s)

  }
)
