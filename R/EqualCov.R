#' covTests
#'
#' Package Documentation
#'
#' @useDynLib covTests
#' @importFrom Rcpp sourceCpp
NULL


.onUnload <- function(libpath){
  library.dynam.unload(
    "covTestR", libpath
  )
}
