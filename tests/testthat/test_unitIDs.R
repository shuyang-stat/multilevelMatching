
context("unitIDs")
## Issue #16


## Matrices

X <- matrix(c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9), ncol=1)
Y <- matrix(c(102,105,120,130,100,80,94,108,96), ncol=1)
# W <- matrix(c(1, 1, 1, 3, 2, 3, 2, 1, 2), ncol=1)
W <- matrix(LETTERS[c(1, 1, 1, 3, 2, 3, 2, 1, 2)], ncol=1)
orig_ids <- paste0("id_", letters[1:NROW(W)])



test_that(
  "Pre-sorting may control the reordering", {

    row.names(W) <- orig_ids
    row.names(Y) <- orig_ids
    row.names(X) <- orig_ids



    # data_list <- list(Y=Y, W=W,X=X)
    data_list <- list(
      Y=Y,
      W=W,
      X=X
    )

    data_ids <- do.call(determineIDs, data_list)
    data_list$unit_ids_unsorted <- orig_ids
    expect_equal( data_list[1:3], data_ids[1:3] )
    expect_equal( data_list$unit_ids_unsorted, data_ids$unit_ids )



    order_W <- order(W[,1])
    data_list <- list(
      Y=Y[order_W,,drop=FALSE],
      W=W[order_W,,drop=FALSE],
      X=X[order_W,,drop=FALSE],
      unit_ids_unsorted = orig_ids[order_W]
    )

    ordered_data <- do.call(reorderByTreatment, data_list)

    expect_equal(
      data_list$X,
      ordered_data$X
    )

    new_W <- as.vector(data_list$W)
    names(new_W) <- rownames(data_list$W)
    expect_equal(
      new_W,
      ordered_data$W
    )

    new_Y <- as.vector(data_list$Y)
    names(new_Y) <- rownames(data_list$Y)
    expect_equal(
      new_Y,
      ordered_data$Y
    )


    # data_list <- list(Y=Y, W=W,X=X)
    data_list <- list(
      Y=Y[order_W,,drop=FALSE],
      W=W[order_W,,drop=FALSE],
      X=X[order_W,,drop=FALSE]
    )

    data_ids <- do.call(determineIDs, data_list)
    data_list$unit_ids_unsorted <- orig_ids[order_W]
    # expect_equal(
    #   data_list,
    #   data_ids
    # )
    expect_equal( data_list[1:3], data_ids[1:3] )
    expect_equal( data_list$unit_ids_unsorted, data_ids$unit_ids )

    data_list <- list(
      Y=Y,
      W=W,
      X=X
    )
    prep_data_args <- append(
      data_list, list( trimming=0, match_on = "multinom", model_options=NULL))
    prep_data_args$M_matches <- 1
    prep_data_args$J_var_matches <- 1
    prep_data <- do.call(prepareData,prep_data_args)

    new_ids <- orig_ids[prep_data$orig_to_sorted]
    expect_equal( new_ids, prep_data$unit_ids_unsorted[prep_data$orig_to_sorted])
    expect_equal( new_ids, prep_data$unit_ids_sorted)
    expect_equal( new_ids, names(prep_data$W) )
    expect_equal( new_ids, names(prep_data$Y) )
    expect_equal( new_ids, rownames(prep_data$X) )

  }
)

test_that(
  "unit IDs are in correct order when W,Y,X are matrices",
  {

    row.names(W) <- orig_ids
    row.names(Y) <- orig_ids
    row.names(X) <- orig_ids
    data_list <- list(Y=Y, W=W,X=X)

    data_ids <- do.call(determineIDs, data_list)
    data_list$unit_ids_unsorted <- orig_ids
    # expect_equal(data_list, data_ids)
    expect_equal( data_list[1:3], data_ids[1:3] )
    expect_equal( data_list$unit_ids_unsorted, data_ids$unit_ids )

    order_W <- order(W[,1])
    data_list <- list(
      Y=Y[order_W,,drop=FALSE],
      W=W[order_W,,drop=FALSE],
      X=X[order_W,,drop=FALSE]
    )
    prep_data_args <- append(
      data_list, list( trimming=0, match_on = "multinom", model_options=NULL))
    prep_data_args$M_matches <- 1
    prep_data_args$J_var_matches <- 1
    prep_data <- do.call(prepareData,prep_data_args)

    new_ids <- orig_ids[order_W][prep_data$orig_to_sorted]
    expect_equal( new_ids, prep_data$unit_ids_unsorted )
    expect_equal( new_ids, names(prep_data$W) )
    expect_equal( new_ids, names(prep_data$Y) )
    expect_equal( new_ids, rownames(prep_data$X) )



    mnom_new <- multiMatch(Y, W, X, trimming=0, match_on = "multinom")

    expect_equal(
      mnom_new$estimate_args$unit_ids_unsorted[mnom_new$estimate_args$orig_to_sorted],
      mnom_new$estimate_args$unit_ids_sorted
    )
    expect_equal(
      mnom_new$estimate_args$unit_ids_sorted[mnom_new$estimate_args$sorted_to_orig],
      mnom_new$estimate_args$unit_ids_unsorted
    )
    expect_equal(
      row.names( mnom_new$impute_mat ),
      mnom_new$estimate_args$unit_ids_unsorted
    )
    expect_equal(
      row.names( mnom_new$impute_mat_sorted ),
      mnom_new$estimate_args$unit_ids_sorted
    )

    expect_equal(
      row.names( mnom_new$impute_mat[mnom_new$estimate_args$orig_to_sorted,] ),
      mnom_new$estimate_args$unit_ids_sorted
    )

    new_ids_sorted <- mnom_new$estimate_args$unit_ids_sorted
    expect_equal(names(mnom_new$estimate_args$Y),new_ids_sorted)
    expect_equal(names(mnom_new$estimate_args$W),new_ids_sorted)
    expect_equal(rownames(mnom_new$estimate_args$X),new_ids_sorted)
    expect_equal(rownames(mnom_new$propensity_scores),new_ids_sorted)



  }
)



test_that(
  "determineIDs() will stop when given duplicates",
  {

    row.names(W) <- orig_ids
    row.names(W)[2] <- orig_ids[1]
    # any(duplicated(row.names(W2)))

    expect_error( determineIDs(Y=Y, W=W, X=X) )
    expect_error( multiMatch(Y, W, X, trimming=0, match_on = "multinom") )

    X <- as.vector(X)
    Y <- as.vector(Y)
    W <- as.vector(W)
    names(Y) <- orig_ids
    names(Y)[2] <- orig_ids[1]

    expect_error( determineIDs(Y=Y, W=W, X=X) )
    expect_error( multiMatch(Y, W, X, trimming=0, match_on = "multinom") )



    names(X) <- orig_ids ## no duplicates in X but there are in Y

    expect_error( determineIDs(Y=Y, W=W, X=X) )
    expect_error( multiMatch(Y, W, X, trimming=0, match_on = "multinom") )

  }
)


test_that(
  "Defensive programming when unit IDs are at odds",
  {

    row.names(W) <- orig_ids
    row.names(Y) <- rev(orig_ids)
    # row.names(X) <- unit_id
    expect_error(determineIDs(Y, W, X))
    expect_error(multiMatch(Y, W, X, trimming=0, match_on = "multinom"))

    X <- as.vector(X)
    expect_error(determineIDs(Y, W, X))
    expect_error(multiMatch(Y, W, X, trimming=0, match_on = "multinom"))

    Y <- as.vector(Y)
    names(X) <- rev(orig_ids)
    expect_error(determineIDs(Y, W, X))
    expect_error(multiMatch(Y, W, X, trimming=0, match_on = "multinom"))

  }
)




test_that(
  "unit IDs are in correct order when only some supplied, when W,Y,X are matrices",{


    row.names(W) <- orig_ids
    row.names(Y) <- NULL
    row.names(X) <- NULL
    data_list <- list(Y=Y, W=W,X=X)

    data_ids <- do.call(determineIDs, data_list)
    data_list$unit_ids_unsorted <- orig_ids

    rownames(data_list$Y) <- orig_ids
    rownames(data_list$X) <- orig_ids
    expect_equal( data_list[1:3], data_ids[1:3] )
    expect_equal( data_list$unit_ids_unsorted, data_ids$unit_ids )

    order_W <- order(W[,1])
    data_list <- list(
      Y=Y[order_W,,drop=FALSE],
      W=W[order_W,,drop=FALSE],
      X=X[order_W,,drop=FALSE]
    )
    prep_data_args <- append(
      data_list, list( trimming=0, match_on = "multinom", model_options=NULL))
    prep_data_args$M_matches <- 1
    prep_data_args$J_var_matches <- 1
    prep_data <- do.call(prepareData,prep_data_args)

    new_ids <- orig_ids[order_W][prep_data$orig_to_sorted]
    expect_equal( new_ids, prep_data$unit_ids_unsorted )
    expect_equal( new_ids, names(prep_data$W) )
    expect_equal( new_ids, names(prep_data$Y) )
    expect_equal( new_ids, rownames(prep_data$X) )



    mnom_new <- multiMatch(Y, W, X, trimming=0, match_on = "multinom")

    expect_equal(
      mnom_new$estimate_args$unit_ids_unsorted[mnom_new$estimate_args$orig_to_sorted],
      mnom_new$estimate_args$unit_ids_sorted
    )
    expect_equal(
      mnom_new$estimate_args$unit_ids_sorted[mnom_new$estimate_args$sorted_to_orig],
      mnom_new$estimate_args$unit_ids_unsorted
    )
    expect_equal(
      row.names( mnom_new$impute_mat ),
      mnom_new$estimate_args$unit_ids_unsorted
    )
    expect_equal(
      row.names( mnom_new$impute_mat_sorted ),
      mnom_new$estimate_args$unit_ids_sorted
    )

    expect_equal(
      row.names( mnom_new$impute_mat[mnom_new$estimate_args$orig_to_sorted,] ),
      mnom_new$estimate_args$unit_ids_sorted
    )

    new_ids_sorted <- mnom_new$estimate_args$unit_ids_sorted
    expect_equal(names(mnom_new$estimate_args$Y),new_ids_sorted)
    expect_equal(names(mnom_new$estimate_args$W),new_ids_sorted)
    expect_equal(rownames(mnom_new$estimate_args$X),new_ids_sorted)
    expect_equal(rownames(mnom_new$propensity_scores),new_ids_sorted)



  }
)

test_that(
  "unit IDs are in correct order when only some supplied, when W and Y are vectors",{

    W <- as.vector(W)
    names(W) <- NULL
    Y <- as.vector(Y)
    names(Y) <- orig_ids
    # row.names(W) <- orig_ids
    # row.names(Y) <- NULL
    row.names(X) <- NULL
    data_list <- list(Y=Y, W=W,X=X)

    data_ids <- do.call(determineIDs, data_list)
    data_list$unit_ids_unsorted <- orig_ids
    names(data_list$W) <- orig_ids
    rownames(data_list$X) <- orig_ids
    # expect_equal(data_list, data_ids)
    expect_equal( data_list[1:3], data_ids[1:3] )
    expect_equal( data_list$unit_ids_unsorted, data_ids$unit_ids )

    order_W <- order(W)
    data_list <- list(
      Y=Y[order_W],
      W=W[order_W],
      X=X[order_W,,drop=FALSE]
    )
    prep_data_args <- append(
      data_list, list( trimming=0, match_on = "multinom", model_options=NULL))
    prep_data_args$M_matches <- 1
    prep_data_args$J_var_matches <- 1
    prep_data <- do.call(prepareData,prep_data_args)

    new_ids <- orig_ids[order_W][prep_data$orig_to_sorted]
    expect_equal( new_ids, prep_data$unit_ids_unsorted )
    expect_equal( new_ids, names(prep_data$W) )
    expect_equal( new_ids, names(prep_data$Y) )
    expect_equal( new_ids, rownames(prep_data$X) )



    mnom_new <- multiMatch(Y, W, X, trimming=0, match_on = "multinom")

    expect_equal(
      mnom_new$estimate_args$unit_ids_unsorted[mnom_new$estimate_args$orig_to_sorted],
      mnom_new$estimate_args$unit_ids_sorted
    )
    expect_equal(
      mnom_new$estimate_args$unit_ids_sorted[mnom_new$estimate_args$sorted_to_orig],
      mnom_new$estimate_args$unit_ids_unsorted
    )
    expect_equal(
      row.names( mnom_new$impute_mat ),
      mnom_new$estimate_args$unit_ids_unsorted
    )
    expect_equal(
      row.names( mnom_new$impute_mat_sorted ),
      mnom_new$estimate_args$unit_ids_sorted
    )
    expect_equal(
      colnames( mnom_new$impute_mat ),
      mnom_new$estimate_args$trt_levels
    )

    expect_equal(
      row.names( mnom_new$impute_mat[mnom_new$estimate_args$orig_to_sorted,] ),
      mnom_new$estimate_args$unit_ids_sorted
    )

    new_ids_sorted <- mnom_new$estimate_args$unit_ids_sorted
    expect_equal(names(mnom_new$estimate_args$Y),new_ids_sorted)
    expect_equal(names(mnom_new$estimate_args$W),new_ids_sorted)
    expect_equal(rownames(mnom_new$estimate_args$X),new_ids_sorted)
    expect_equal(rownames(mnom_new$propensity_scores),new_ids_sorted)



  }
)




