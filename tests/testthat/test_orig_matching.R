
context("test_orig matching output from v0.1.0")

X <- c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9)
Y <- c(102,105,120,130,100,80,94,108,96)
W <- c(1,1,1,3,2,3,2,1,2)
names(Y) <- letters[2+(13:5)]

test_that(
  "matchX returns same output as original",  {
    set.seed(22)
    fit <- multilevelMatchX(Y,W,X)

    fit_orig <- structure(list(
      tauestimate = structure(
        c(-10.6666666666667, 6.66666666666666, 17.3333333333333),
        .Names = c("EY(2)-EY(1)", "EY(3)-EY(1)", "EY(3)-EY(2)")),
      varestimate = structure(
        c(9.11111111111111, 615.58024691358, 613.925925925926),
        .Names = c("EY(2)-EY(1)", "EY(3)-EY(1)", "EY(3)-EY(2)"))),
      .Names = c("tauestimate", "varestimate"))


    expect_equal(
      fit$tauestimate,
      fit_orig$tauestimate,
      tol=1e-5
    )

    expect_equal(
      fit$varestimate,
      fit_orig$varestimate,
      tol=1e-5
    )

    if (names(fit)[[1]] == "tauestimate"){
      ## original code from v0.1.0
      expect_equal(
        fit$tauestimate,
        fit_orig$tauestimate,
        tol=1e-5
      )

      expect_equal(
        fit$varestimate,
        fit_orig$varestimate,
        tol=1e-5
      )
    } else {
      expect_equal(
        fit$results$Estimate,
        fit_orig$tauestimate,
        tol=1e-5,
        check.names = FALSE
      )

      expect_equal(
        fit$results$Variance,
        fit_orig$varestimate,
        tol=1e-5,
        check.names = FALSE
      )
    }

  }
)
test_that(
  "GPSMatch on MLR returns same output as original",  {
    set.seed(22)
    fit <- multilevelGPSMatch(Y,W,X,Trimming=0,GPSM="multinomiallogisticReg")
    set.seed(22)
    fit2 <- multilevelGPSMatch(Y,W,X,Trimming=1,GPSM="multinomiallogisticReg")


    fit_orig <-
      structure( list(
          tauestimate = structure(
            c(-10.4444444444444, 6.66666666666666,17.1111111111111),
            .Names = c("EY(2)-EY(1)", "EY(3)-EY(1)", "EY(3)-EY(2)")
          ),
          varestimate = structure(
            c(8.54595336076818, 616.913580246914, 611.122085048011),
            .Names = c("EY(2)-EY(1)", "EY(3)-EY(1)", "EY(3)-EY(2)")
          ),
          varestimateAI2012 = structure(
            c(8.30202363928055, 411.456234042855, 434.24703718368),
            .Names = c("EY(2)-EY(1)", "EY(3)-EY(1)", "EY(3)-EY(2)")
          ),analysisidx = 1:9
        ),.Names = c("tauestimate","varestimate","varestimateAI2012","analysisidx"
        )
      )

    fit2_orig <-

      structure(
        list(
          tauestimate = structure(
            c(-9.375, 5.875, 15.25),
            .Names = c("EY(2)-EY(1)", "EY(3)-EY(1)", "EY(3)-EY(2)")
          ),
          varestimate = structure(
            c(7.794921875, 582.654296875, 576.3046875),
            .Names = c("EY(2)-EY(1)", "EY(3)-EY(1)", "EY(3)-EY(2)")
          ),
          varestimateAI2012 = structure(
            c(5.07205734159109, 383.848574889225,

              430.978089158589),
            .Names = c("EY(2)-EY(1)", "EY(3)-EY(1)",

                       "EY(3)-EY(2)")
          ),
          analysisidx = structure(
            c(1L, 2L, 4L, 5L,
              6L, 7L, 8L, 9L),
            .Names = c("1", "2", "4", "5", "6", "7",
                       "8", "9")
          )
        ),
        .Names = c("tauestimate",
                   "varestimate",
                   "varestimateAI2012",
                   "analysisidx"))



    expect_equal(
      fit$tauestimate,
      fit_orig$tauestimate,
      tol=1e-5
    )

    expect_equal(
      fit$varestimate,
      fit_orig$varestimate,
      tol=1e-5
    )
    expect_equal(
      fit$varestimateAI2012,
      fit_orig$varestimateAI2012,
      tol=1e-5
    )

    expect_equal(
      fit$analysisidx,
      fit_orig$analysisidx,
      tol=1e-5
    )


    expect_equal(
      fit2$tauestimate,
      fit2_orig$tauestimate,
      tol=1e-5
    )

    expect_equal(
      fit2$varestimate,
      fit2_orig$varestimate,
      tol=1e-5
    )
    expect_equal(
      fit2$varestimateAI2012,
      fit2_orig$varestimateAI2012,
      tol=1e-5
    )

    expect_equal(
      fit2$analysisidx,
      fit2_orig$analysisidx,
      tol=1e-5
    )

    if (names(fit)[[1]] == "tauestimate"){
      ## original code from v0.1.0
      expect_equal(
        fit$tauestimate,
        fit_orig$tauestimate,
        tol=1e-5
      )

      expect_equal(
        fit$varestimate,
        fit_orig$varestimate,
        tol=1e-5
      )
      expect_equal(
        fit$varestimateAI2012,
        fit_orig$varestimateAI2012,
        tol=1e-5
      )

      expect_equal(
        fit$analysisidx,
        fit_orig$analysisidx,
        tol=1e-5
      )


      expect_equal(
        fit2$tauestimate,
        fit2_orig$tauestimate,
        tol=1e-5
      )

      expect_equal(
        fit2$varestimate,
        fit2_orig$varestimate,
        tol=1e-5
      )
      expect_equal(
        fit2$varestimateAI2012,
        fit2_orig$varestimateAI2012,
        tol=1e-5
      )

      expect_equal(
        fit2$analysisidx,
        fit2_orig$analysisidx,
        tol=1e-5
      )


    } else {
      expect_equal(
        fit$results$Estimate,
        check.names = FALSE,
        fit_orig$tauestimate,
        tol=1e-5
      )

      expect_equal(
        fit$results$Variance,
        check.names = FALSE,
        fit_orig$varestimate,
        tol=1e-5
      )
      expect_equal(
        fit$results$VarianceAI2012,
        check.names = FALSE,
        fit_orig$varestimateAI2012,
        tol=1e-5
      )

      ## new code outputs NULL for analysis_idx when no trimming
      expect_failure(
        expect_equal(
          fit$analysis_idx$indices_kept,
          # check.names = FALSE,
          fit_orig$analysisidx,
          tol=1e-5
        )
      )


      expect_equal(
        fit2$results$Estimate,
        check.names = FALSE,
        fit2_orig$tauestimate,
        tol=1e-5
      )

      expect_equal(
        fit2$results$Variance,
        check.names = FALSE,
        fit2_orig$varestimate,
        tol=1e-5
      )
      expect_equal(
        fit2$results$VarianceAI2012,
        check.names = FALSE,
        fit2_orig$varestimateAI2012,
        tol=1e-5
      )

      expect_equal(
        fit2$analysis_idx$indices_kept,
        fit2_orig$analysisidx,
        tol=1e-5
      )
    }


  }
)

test_that(
  "GPSMatch on POLR returns same output as original",  {
    set.seed(22)
    fit <- multilevelGPSMatch(Y,W,X,Trimming=0,GPSM="ordinallogisticReg")
    set.seed(22)
    fit2 <- multilevelGPSMatch(Y,W,X,Trimming=1,GPSM="ordinallogisticReg")


    fit_orig <-structure(list(tauestimate = structure(c(-10.6666666666667, 6.66666666666666,
                                                        17.3333333333333), .Names = c("EY(2)-EY(1)", "EY(3)-EY(1)", "EY(3)-EY(2)"
                                                        )), varestimate = structure(c(9.11111111111111, 615.58024691358,
                                                                                      613.925925925926), .Names = c("EY(2)-EY(1)", "EY(3)-EY(1)", "EY(3)-EY(2)"
                                                                                      )), varestimateAI2012 = structure(c(NA, NA, NA), .Names = c("EY(2)-EY(1)",
                                                                                                                                                  "EY(3)-EY(1)", "EY(3)-EY(2)")), analysisidx = 1:9), .Names = c("tauestimate",
                                                                                                                                                                                                                 "varestimate", "varestimateAI2012", "analysisidx"))

    fit2_orig <-
      structure(
        list(
          tauestimate = structure(
            c(-9, 6.25, 15.25),
            .Names = c("EY(2)-EY(1)",
                       "EY(3)-EY(1)", "EY(3)-EY(2)")
          ),
          varestimate = structure(
            c(7.90625,
              583.5859375, 576.3046875),
            .Names = c("EY(2)-EY(1)", "EY(3)-EY(1)",
                       "EY(3)-EY(2)")
          ),
          varestimateAI2012 = structure(
            c(NA, NA, NA),
            .Names = c("EY(2)-EY(1)",
                       "EY(3)-EY(1)", "EY(3)-EY(2)")
          ),
          analysisidx = structure(
            c(1L,
              2L, 4L, 5L, 6L, 7L, 8L, 9L),
            .Names = c("1", "2", "4", "5", "6",
                       "7", "8", "9")
          )
        ),
        .Names = c("tauestimate", "varestimate", "varestimateAI2012",
                   "analysisidx"))


    if (names(fit)[[1]] == "tauestimate"){

    expect_equal(
      fit$tauestimate,
      fit_orig$tauestimate,
      tol=1e-5
    )

    expect_equal(
      fit$varestimate,
      fit_orig$varestimate,
      tol=1e-5
    )
    expect_equal(
      fit$varestimateAI2012,
      fit_orig$varestimateAI2012,
      tol=1e-5
    )

    expect_equal(
      fit$analysisidx,
      fit_orig$analysisidx,
      tol=1e-5
    )


    expect_equal(
      fit2$tauestimate,
      fit2_orig$tauestimate,
      tol=1e-5
    )

    expect_equal(
      fit2$varestimate,
      fit2_orig$varestimate,
      tol=1e-5
    )
    expect_equal(
      fit2$varestimateAI2012,
      fit2_orig$varestimateAI2012,
      tol=1e-5
    )

    expect_equal(
      fit2$analysisidx,
      fit2_orig$analysisidx,
      tol=1e-5
    )


    } else {
      expect_equal(
        fit$results$Estimate,
        check.names = FALSE,
        fit_orig$tauestimate,
        tol=1e-5
      )

      expect_equal(
        fit$results$Variance,
        check.names = FALSE,
        fit_orig$varestimate,
        tol=1e-5
      )
      expect_equal(
        fit$results$VarianceAI2012,
        check.names = FALSE,
        fit_orig$varestimateAI2012,
        tol=1e-5
      )

      ## new code outputs NULL for analysis_idx when no trimming
      expect_failure(
        expect_equal(
          fit$analysis_idx$indices_kept,
          # check.names = FALSE,
          fit_orig$analysisidx,
          tol=1e-5
        )
      )


      expect_equal(
        fit2$results$Estimate,
        check.names = FALSE,
        fit2_orig$tauestimate,
        tol=1e-5
      )

      expect_equal(
        fit2$results$Variance,
        check.names = FALSE,
        fit2_orig$varestimate,
        tol=1e-5
      )
      expect_equal(
        fit2$results$VarianceAI2012,
        check.names = FALSE,
        fit2_orig$varestimateAI2012,
        tol=1e-5
      )

      expect_equal(
        fit2$analysis_idx$indices_kept,
        fit2_orig$analysisidx,
        tol=1e-5
      )
    }


  }
)
