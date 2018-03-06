
context("calcKMVarFactor")

test_that(
  "calcKMVarFactor returns the correct values of K_M_var_factor from Kiw",{

    ## Expected vectors were calculated by hand from the components presented in
    ## Theorem 7 in Abadie and Imbens 2006 Econometrica paper

    expect_equal(
      calcKMVarFactor(Kiw = c(1,5,10,23), M_matches = 1),
      c(2, 30, 110, 552)
    )
    expect_equal(
      calcKMVarFactor(Kiw = c(1,5,10,23), M_matches = 2),
      c(1, 10, 32.5, 149.5)
    )
}
)
