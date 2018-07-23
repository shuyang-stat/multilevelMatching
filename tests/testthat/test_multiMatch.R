
context("test_multiMatch() versus legacy matching funs")


X <- matrix(c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9), ncol=1)
Y <- matrix(c(102,105,120,130,100,80,94,108,96), ncol=1)
W <- matrix(c(1, 1, 1, 3, 2, 3, 2, 1, 2), ncol=1)


existing_GPS_matrix <- cbind(
  c(0.5, 0.3, 0.5, 0.5, 0.5, 0.3, 0.3, 0.5, 0.3),
  c(1,1.6, 1, 1, 1, 1.6,1.6, 1,1.6)/6,
  c(2, 2.6, 2, 2, 2, 2.6, 2.6, 2, 2.6)/6
)
eps_legacy <- multilevelGPSMatch(
  Y = Y,
  W = W,
  X = existing_GPS_matrix,
  Trimming = 0,
  GPSM = "existing"
)

expect_message(
  eps_multiMatch <- multiMatch(
    Y = Y,
    W = W,
    X = existing_GPS_matrix,
    trimming = 0,
    match_on = "existing"
  )
)


##  Matching on "existing" is not passing when ties are present
##  For more tests, see test_existing_GPS_matching.R script.
test_that(
  "Discrepancy when matching on existing GPS between multiMatch() and multilvelGPSMatch()",
  {


    expect_identical(
      object = names(eps_legacy$varestimate),
      expected = (eps_multiMatch$results)$Param
    )
    expect_failure(
      expect_equal(
        object = eps_legacy$tauestimate,
        expected = (eps_multiMatch$results)$Estimate,
        tolerance = 1e-7,
        check.names = FALSE
      )
    )
    expect_failure(
      expect_equal(
        object = eps_legacy$varestimate,
        expected = (eps_multiMatch$results)$Variance,
        tolerance = 1e-7,
        check.names = FALSE
      )
    )

  }
)


## not passing tests
test_that("match on GPS with existing GPS DOES NOT RETURN same output", {

  # Fails testthat::test() but passes devtools::check() ??
  # expect_failure(
  #   expect_equal( (eps_legacy$results)$Estimate,
  #                 c( -8.000000 , 1.777778 , 9.777778),
  #                 tolerance = my_tolerance)
  # )

  expect_failure(
    expect_equal( (eps_legacy$results)$Variance,
                  c( 18.04938 ,573.25377, 552.78464),
                  tolerance = my_tolerance)
  )
})

test_that(
  "multiMatch() on existing GPS has discrepant results when ties present",
  {
    expect_failure(
      expect_equal(
        (eps_multiMatch$results)$Estimate,
        c( -8.000000 , 1.777778 , 9.777778),
        tolerance = my_tolerance
      )
    )

    expect_failure(
      expect_equal(
        (eps_multiMatch$results)$Variance,
        c( 18.04938 ,573.25377, 552.78464),
        tolerance = my_tolerance
      )
    )
  }
)



test_that(
  "multinom-matching works same between multiMatch() and legacy",
  {

    ## old function
    mnom_legacy <- multilevelGPSMatch(Y,W,X,Trimming=0,GPSM="multinomiallogisticReg")

    # new function
    expect_message(
      mnom_multiMatch <-
        multiMatch(Y, W, X, trimming=0, match_on = "multinom")
    )

    expect_identical(
      object = names(mnom_legacy$varestimate),
      expected = (mnom_multiMatch$results)$Param
    )

    expect_equal(
      object = mnom_legacy$tauestimate,
      expected = (mnom_multiMatch$results)$Estimate,
      tolerance = 1e-7,
      check.names = FALSE
    )

    expect_equal(
      object = mnom_legacy$varestimate,
      expected = (mnom_multiMatch$results)$Variance,
      tolerance = 1e-7,
      check.names = FALSE
    )
    expect_equal(
      object = mnom_legacy$varestimateAI2012,
      expected = (mnom_multiMatch$results)$VarianceAI2016,
      tolerance = 1e-7,
      check.names = FALSE
    )
  }
)






test_that(
  "polr-matching works same between multiMatch() and legacy",
  {

    ##old
    polr_legacy <- multilevelGPSMatch(Y,W,X,Trimming=0,GPSM="ordinallogisticReg")
    #new
    expect_message(
      polr_multiMatch <-
        multiMatch(Y, W, X, trimming=0, match_on = "polr")
    )

    expect_identical(
      object = names(polr_legacy$varestimate),
      expected = (polr_multiMatch$results)$Param
    )

    expect_equal(
      object = polr_legacy$tauestimate,
      expected = (polr_multiMatch$results)$Estimate,
      tolerance = 1e-7,
      check.names = FALSE
    )

    expect_equal(
      object = polr_legacy$varestimate,
      expected = (polr_multiMatch$results)$Variance,
      tolerance = 1e-7,
      check.names = FALSE
    )
    expect_equal(
      object = polr_legacy$varestimateAI2012,
      expected = (polr_multiMatch$results)$VarianceAI2016,
      tolerance = 1e-7,
      check.names = FALSE
    )
  }
)





test_that(
  "covariate-matching works same between multiMatch() and legacy",
  {
    ## old
    covar_legacy <- multilevelMatchX(Y=Y,W=W,X=X)

    ## new
    expect_message(
      covar_multiMatch <-
        multiMatch(Y=Y,W=W,X=X,trimming=0,match_on="covariates")
    )

    expect_identical(
      object = names(covar_legacy$varestimate),
      expected = (covar_multiMatch$results)$Param
    )

    expect_equal(
      object = covar_legacy$tauestimate,
      expected = (covar_multiMatch$results)$Estimate,
      tolerance = 1e-7,
      check.names = FALSE
    )

    expect_equal(
      object = covar_legacy$varestimate,
      expected = (covar_multiMatch$results)$Variance,
      tolerance = 1e-7,
      check.names = FALSE
    )
  }
)

