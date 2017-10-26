#' covTests
#'
#' Package Documentation
#'
#' @useDynLib covTestR
#' @importFrom Rcpp sourceCpp
NULL


.onUnload <- function(libpath){
  library.dynam.unload(
    "covTestR", libpath
  )
}

