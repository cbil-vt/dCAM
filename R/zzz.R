#' @importFrom rJava .jpackage .jclassPath
.onLoad <- function(libname, pkgname) {
  .jpackage(pkgname, lib.loc = libname)
  cat("my Java class path: ")
  print(.jclassPath())
  # print(.jclassPath(.jclassLoader(package=pkgname)))
}
