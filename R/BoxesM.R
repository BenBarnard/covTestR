#' Test of Equality of Covariances given by Box's M
#'
#' @inheritParams Chaipitak2013
#'
#' @return Test statistic for Chaipitak 2013
#'
#' @export
#'
#' @examples BoxsM(iris, group = Species)
#'
BoxesM <- function(x, ...){
  ls <- lazy_dots(...)
  matrix_ls <- x

  statistic <- BoxesMStat(matrix_ls)


  out <- list(statistic = statistic)
  out
}

