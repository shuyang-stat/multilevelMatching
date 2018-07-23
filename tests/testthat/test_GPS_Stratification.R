
context("Stratify on GPS")

load(   # Y,W,X,  match4,  match4_lp,
  file = quickLookup("test_stratification_data.Rdata")
)

test_that(
  "GPS stratification returns original (v0.1) results",
  {
    set.seed(22)
    strat <- multilevelGPSStratification(
      Y,W,X,NS=10,GPSM="multinomiallogisticReg",
      linearp=0,nboot=5
    )

    expect_equal(
      strat$tauestimate,
      match4$tauestimate,
      tolerance = my_tolerance,
      check.attributes = FALSE
    )
    expect_equal(
      strat$varestimate,
      match4$varestimate,
      tolerance = my_tolerance,
      check.attributes = FALSE
    )
    expect_identical(
      names(strat$tauestimate),
      names(match4$varestimate)
    )


    set.seed(22)
    strat_lp <- multilevelGPSStratification(
      Y,W,X,NS=10,
      GPSM="ordinallogisticReg",linearp=1,nboot=5
    )

    expect_equal(
      strat_lp$tauestimate,
      match4_lp$tauestimate,
      tolerance = my_tolerance,
      check.attributes = FALSE
    )
    expect_equal(
      strat_lp$varestimate,
      match4_lp$varestimate,
      tolerance = my_tolerance,
      check.attributes = FALSE
    )
    expect_identical(
      names(strat_lp$tauestimate),
      names(match4_lp$varestimate)
    )
  }
)
