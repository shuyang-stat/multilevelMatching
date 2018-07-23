
context("test_multiple_matches: M, J >=1")


test_that(
  "Imputed outcomes same as before for M=1",{

    load(file = quickLookup("test_M_match.Rdata"))

    #                                    #
    # Step 3: Calculate some components  #
    #                                    #
    set.seed(22)
    diff_trt_matching_args$M <- 1
    diff_trt_match <- do.call(Matching::Match, diff_trt_matching_args)



    est_kk_list <- wrangleImputations(
      match_output = diff_trt_match,
      M_matches = diff_trt_matching_args$M,
      Y = Y,
      which_same_trt = 1:4,
      which_diff_trt = 5:9
    )


    expect_equal(
      est_kk_list$Yiw_kk,
      c(102, 105, 120, 108, 102, 105, 108, 105, 105)
    )

    expect_equal(
      est_kk_list$mean_po_kk,
      106.666667,
      tol=1e-4
    )

    expect_equal(
      as.numeric(est_kk_list$Kiw),
      c(1,3,0,1),
      check.attributes = FALSE
    )

  }
)

test_that(
  "Imputed outcomes same as expected for M=2",{


    load(file = quickLookup("test_M_match.Rdata"))

    #                                    #
    # Step 3: Calculate some components  #
    #                                    #
    set.seed(22)
    diff_trt_matching_args$M <- 2
    diff_trt_match <- do.call(Matching::Match, diff_trt_matching_args)


    est_kk_list <- wrangleImputations(
      match_output = diff_trt_match,
      M_matches = diff_trt_matching_args$M,
      Y = Y,
      which_same_trt = 1:4,
      which_diff_trt = 5:9
    )

    expect_equal(
      est_kk_list$Yiw_kk,
      c(102, 105, 120, 108, 105, 103.5, 105, 103.5, 103.5)
    )

    expect_equal(
      est_kk_list$mean_po_kk,
      106.16666667,
      tol=1e-4
    )

    expect_equal(
      as.numeric(est_kk_list$Kiw),
      c(5,3,0,2),
      check.attributes = FALSE
    )
  }
)

## TODO:: test that this code works when Matching drops units


test_that(
  "multiMatch() does not throw errors for M_matches>=1",{


    ### adding three more individuals so that J=2 doesn't throw errors
    X <- matrix(c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9,10,10,10), ncol=1)
    Y <- matrix(c(102,105,120,130,100,80,94,108,96,100,100,100), ncol=1)
    W <- matrix(c(1, 1, 1, 3, 2, 3, 2, 1, 2,1,2,3), ncol=1)
    rownames(W) <- letters[1:12]

    expect_silent(
      output <- multiMatch(
        Y, W, X,
        trimming = 0,
        match_on = "multinom",
        M_matches = 1
      )
    )


    expect_silent(
      output <- multiMatch(
        Y, W, X,
        trimming = 0,
        match_on = "multinom",
        M_matches = 2
      )
    )

    expect_silent(
      output <- multiMatch(
        Y, W, X,
        trimming = 0,
        match_on = "multinom",
        M_matches = 2,
        J_var_matches = 2
      )
    )


  }
)



test_that(
  "calcSigSqAI2006() returns same output as before for J=1",{

    same_trt_match <- readRDS(file =  quickLookup("J_var_matches_1_match.Rds"))

    same_trt_match_data <- same_trt_match$mdata
    J_var_matches <- 1
    sigsqiw_old <- ((J_var_matches)/(1+J_var_matches))*
      (same_trt_match_data$Y[which(same_trt_match_data$Tr==1)] -
         same_trt_match_data$Y[which(same_trt_match_data$Tr==0)])^2

    sigsqiw_new <- calcSigSqAI2006(
      match_output  = same_trt_match, J = J_var_matches
    )

    expect_equal(sigsqiw_old,sigsqiw_new)
  }
)


test_that(
  "calcSigSqAI2006() returns correct output for J=2",{

    same_trt_match <- readRDS(file =  quickLookup("J_var_matches_2_match.Rds"))
    same_trt_match_data <- same_trt_match$mdata
    J_var_matches <- 2

    temp_orig_outcome <- same_trt_match_data$Y[which(same_trt_match_data$Tr==1)]
    orig_outcomes <- temp_orig_outcome[c(1,3,5,7)] ## theyre repeated

    temp_matched_outcomes <- same_trt_match_data$Y[which(same_trt_match_data$Tr==0)]
    n_here <- length(orig_outcomes)
    matched_outcomes_averaged <- rep(NA, n_here)
    for (ii in 1:n_here){
      indices_ii <- (2*ii) + (-1:0)
      matched_outcomes_averaged[ii] <- mean(
        temp_matched_outcomes[indices_ii]
      )
    }

    sigsqiw_by_hand <- ((J_var_matches)/(1+J_var_matches)) *
      ( orig_outcomes - matched_outcomes_averaged )^2

    sigsqiw_new <- calcSigSqAI2006(
      match_output  = same_trt_match, J = J_var_matches
    )

    expect_equal(sigsqiw_by_hand,sigsqiw_new)
  }
)

test_that(
  "multiMatch() does not throw errors for J_var_matches>=1",{


    ### adding three more individuals so that J=2 doesn't throw errors
    X <- matrix(c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9,10,10,10), ncol=1)
    Y <- matrix(c(102,105,120,130,100,80,94,108,96,100,100,100), ncol=1)
    W <- matrix(c(1, 1, 1, 3, 2, 3, 2, 1, 2,1,2,3), ncol=1)
    rownames(W) <- letters[1:12]

    expect_silent(
      output <- multiMatch(
        Y, W, X,
        trimming=0,
        match_on = "multinom",
        J_var_matches=1
      )
    )


    expect_silent(
      output <- multiMatch(
        Y, W, X,
        trimming=0,
        match_on = "multinom",
        J_var_matches=2
      )
    )

  }
)

test_that(
  "catching errors for bad values of J_var_matches",{

    X <- matrix(c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9), ncol=1)
    Y <- matrix(c(102,105,120,130,100,80,94,108,96), ncol=1)
    W <- matrix(c(1, 1, 1, 3, 2, 3, 2, 1, 2), ncol=1)
    rownames(W) <- letters[1:9]
        ## Not enough inidividuals to match to
    ##    there are only 2 individuals with W=3,
    ##    so you can't find J=2 matches for them.
    ##    Matching::Match would throw a warning,
    ##    so multiMatch() should throw an error.
    expect_error(
      # expect_warning(
       output <- multiMatch(
        Y, W, X,
        trimming=0,
        match_on = "multinom",
        J_var_matches=2
       )
      # )
    )

    expect_error(
      # expect_warning(
       output <- multiMatch(
        Y, W, X,
        trimming=0,
        match_on = "multinom",
        M_matches = 5,
        J_var_matches=1
       )
      # )
    )

    ### adding three more individuals so that J=2 doesn't throw errors
    X <- matrix(c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9,10,10,10), ncol=1)
    Y <- matrix(c(102,105,120,130,100,80,94,108,96,100,100,100), ncol=1)
    W <- matrix(c(1, 1, 1, 3, 2, 3, 2, 1, 2,1,2,3), ncol=1)
    rownames(W) <- letters[1:12]

    ## Values of M that aren't supported
    expect_error(
       output <- multiMatch(
        Y, W, X,
        trimming=0,
        match_on = "multinom",
        M_matches=1.5 ## not a matchable number
       )
    )

    expect_error(
       output <- multiMatch(
        Y, W, X,
        trimming=0,
        match_on = "multinom",
        M_matches=0 ## not a matchable number
       )
    )
    expect_error(
       output <- multiMatch(
        Y, W, X,
        trimming=0,
        match_on = "multinom",
        M_matches=-2 ## not a matchable number
       )
    )
    ## Values of J that aren't supported
    expect_error(
       output <- multiMatch(
        Y, W, X,
        trimming=0,
        match_on = "multinom",
        J_var_matches=1.5 ## not a matchable number
       )
    )

    expect_error(
       output <- multiMatch(
        Y, W, X,
        trimming=0,
        match_on = "multinom",
        J_var_matches=0 ## not a matchable number
       )
    )
    expect_error(
       output <- multiMatch(
        Y, W, X,
        trimming=0,
        match_on = "multinom",
        J_var_matches=-2 ## not a matchable number
       )
    )

    }
)

