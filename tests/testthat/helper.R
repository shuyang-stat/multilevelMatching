

quickLookup <- function(name) {
  rprojroot::find_testthat_root_file("testing_datafiles", name)
}
