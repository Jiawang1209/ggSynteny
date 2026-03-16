.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste(
      "ggSynteny v0.1.0",
      "Genome synteny data processing and visualization",
      sep = "\n"
    )
  )
}
