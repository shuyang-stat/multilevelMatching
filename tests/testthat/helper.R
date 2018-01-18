

quickLookup <- function(name) {
  rprojroot::find_package_root_file(
    "tests", "testthat", "testing_datafiles", name)
}

my_tolerance <- 1e-7
