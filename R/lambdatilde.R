#' Noise reduction lambda
#'
#' @param x
#' @param y
#'
#' @export
#'
lambdatilde <- function(x, y){
  lambdaest <- sapply(1:length(x), function(z){
    x[[z]] - (tr(y) - sum(x)) / (nrow(y) - 1 - z)
  })
}
