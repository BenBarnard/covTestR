#' covTestR
#'
#' Package Documentation
#'
#' @useDynLib covTestR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL


.onUnload <- function(libpath){
  library.dynam.unload(
    "covTestR", libpath
  )
}

