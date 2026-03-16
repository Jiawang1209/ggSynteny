extdata_path <- function(filename) {
  installed_path <- system.file("extdata", filename, package = "ggSynteny")
  if (nzchar(installed_path)) {
    return(installed_path)
  }

  testthat::test_path("../../inst/extdata", filename)
}
