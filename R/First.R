#  First.R
#
#  $Revision: 1.1 $ $Date: 2013/10/19 03:06:59 $
#

.onLoad <- function(...) {} 

.onAttach <- function(libname, pkgname) {
  dfile <- system.file("DESCRIPTION", package="polyclip")
  vs <- read.dcf(file=dfile, fields="Version")
  cl <- read.dcf(file=dfile, fields="ClipperInfo")
  msg <- paste("polyclip", vs, "built from clipper", cl)
  packageStartupMessage(msg)
  invisible(NULL)
}

  
