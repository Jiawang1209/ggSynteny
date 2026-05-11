.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste(
      "ggSynteny v0.1.1",
      "Genome synteny data processing and visualization",
      sep = "\n"
    )
  )
}
