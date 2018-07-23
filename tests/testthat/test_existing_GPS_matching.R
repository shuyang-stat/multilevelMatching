
context("Matching on existing GPS works (with limited discrepancy)")
# skip("skip this test for now")


# X<-c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9)
Y<-c(102,105,120,130,100,80,94,108,96)
W<-c(1,1,1,3,2,3,2,1,2)


existing_GPS_matrix <- cbind(
  c(0.5, 0.3, 0.5, 0.5, 0.5, 0.3, 0.3, 0.5, 0.3),
  c(1,1.6, 1, 1, 1, 1.6,1.6, 1,1.6)/6,
  c(2, 2.6, 2, 2, 2, 2.6, 2.6, 2, 2.6)/6
)
rownames(existing_GPS_matrix) <- LETTERS[6+(10:2)]

set.seed(12345)

run_legacy <- multilevelGPSMatch(
  Y=Y,W=W,X=existing_GPS_matrix,
  Trimming=0,GPSM="existing")

test_that(
  "multilevelGPSMatch (with existing GPS) returns original (v0.1) output",
  {

    expect_equal(
      run_legacy$tauestimate,
      c(-8.888889, 1.111111, 10.000000),
      tolerance = my_tolerance,
      check.names = FALSE
    )
    expect_equal(
      run_legacy$varestimate,
      c(26.01097, 578.20850, 551.50617),
      tolerance = my_tolerance,
      check.names = FALSE
    )
  }
)


set.seed(12345)
run_multiMatch <- multiMatch(
  Y=Y,
  W=W,
  X=existing_GPS_matrix,
  trimming=0,
  match_on="existing"
)

## Some tests fail!
##    I think this is when there are ties to break
##    because the order of the matching process is different in multiMatch
test_that(
  "Discrepancy when matching on existing GPS with ties",
  {

    expect_identical(
      object = names(run_legacy$varestimate),
      expected = (run_multiMatch$results)$Param
    )
    expect_failure(
      expect_equal(
        run_legacy$tauestimate,
        (run_multiMatch$results)$Estimate,
        tolerance = my_tolerance,
        check.attributes = FALSE
      )
    )
    expect_failure(
      expect_equal(
        run_legacy$varestimate,
        (run_multiMatch$results)$Variance,
        tolerance = my_tolerance,
        check.attributes = FALSE
      )
    )
  }
)


## multiMatch has been stable for a while though.
test_that(
  "multiMatch() on GPS with existing GPS returns same values as it used to",
  {
    run_multiMatch_orig <- readRDS(file = quickLookup("existingGPS_t4mm_orig.Rds"))

    names(run_multiMatch_orig$results)[6] <- "VarianceAI2016" #2018-06-10

    expect_equal(
      (run_multiMatch_orig$results) ,
      (run_multiMatch$results) ,
      tolerance = my_tolerance
    )
  }
)


## Existing GPS matching via multiMatch() is not perfect, but close

test_that(
  "multiMatch() on existing GPS returns SIMILAR to multilevelGPSMatch() from v0.1",
  {
    set.seed(11)
    N <- 300
    X <- rnorm(N, 0, 1)
    prW1 <- sample(size = N, x=(1:4)/10, replace = TRUE)
    prW2 <- (1-prW1)*sample(size = N, x=(1:4)/5, replace = TRUE)
    prW3 <- 1- (prW1+prW2)
    existing_GPS_matrix <- cbind(prW1, prW2, prW3)
    W <- rep(NA, N)
    for(ii in 1:N){
      W[ii] <- sample(1:3, size = 1, replace = TRUE,
                      prob = existing_GPS_matrix[ii,])
    }
    Y <- round(rnorm(N, 10 - W +0.2*X, 1),1)

    # existing_GPS_matrix <- cbind(pr_w1, pr_w2,pr_w3)

    set.seed(12345)
    run_legacy <- multilevelGPSMatch(
      Y = Y,
      W = W,
      X = existing_GPS_matrix,
      Trimming = 0,
      GPSM = "existing"
    )

    set.seed(12345)
    expect_message(
      run_multiMatch <- multiMatch(
        Y = Y,
        W = W,
        X = existing_GPS_matrix,
        trimming = 0,
        match_on = "existing"
      )
    )

    expect_identical(
      object = names(run_legacy$varestimate),
      expected = (run_multiMatch$results)$Param
    )
    expect_equal(
      run_legacy$tauestimate,
      (run_multiMatch$results)$Estimate,
      tolerance = 0.04,
      check.names = FALSE
    )
    expect_equal(
      run_legacy$varestimate,
      (run_multiMatch$results)$Variance,
      tolerance = 0.008,
      check.attributes = FALSE
    )
  }
)



## When there are no ties, the multiMatch() may return the same thing as orig

test_that(
  paste(
    "multiMatch on GPS with existing GPS agrees with multilevelGPSMatch",
    "output when there are no ties to break in the existing GPS"
  ), {




    eX <- existing_GPS_matrix+rnorm(27,0,0.01)
    eX <- eX / rowSums(eX)
    set.seed(12345)

    run_multiMatch_eX <- multiMatch(
      Y=Y,W=W,X=eX,
      trimming=0,match_on="existing"
    )


    set.seed(12345)
    run_legacy_eX <- multilevelGPSMatch(
      Y=Y,W=W,X=eX,Trimming=0,GPSM="existing"
    )

    expect_equal(
      run_legacy_eX$tauestimate,
      (run_multiMatch_eX$results)$Estimate,
      tolerance = my_tolerance,
      check.names = FALSE
    )
    expect_equal(
      run_legacy_eX$varestimate,
      (run_multiMatch_eX$results)$Variance,
      tolerance = my_tolerance,
      check.names = FALSE
    )

})
