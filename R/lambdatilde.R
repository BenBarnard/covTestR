#' Noise reduction lambda
#'
#' @param x x
#' @param y y
#'
#' @keywords internal
#'
#'
lambdatilde <- function(x, y){
  lambdaest <- sapply(1:length(x), function(z){
    x[[z]] - (tr(y) - sum(x)) / (nrow(y) - 1 - z)
  })
}
