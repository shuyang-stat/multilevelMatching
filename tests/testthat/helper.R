

quickLookup <- function(name) {
  # rprojroot::find_package_root_file(
  #   "tests", "testthat",
  rprojroot::find_testthat_root_file(
    "testing_datafiles", name)
}

my_tolerance <- 1e-7
