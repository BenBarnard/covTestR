#' Test of Equality of Covariances given by Box's M
#'
#' @inheritParams Chaipitak2013
#'
#' @return Test statistic for Chaipitak 2013
#'
#' @export
#'
#' @examples 
#' irisSpecies <- unique(iris$Species)
#' 
#' iris_ls <- lapply(irisSpecies, 
#'     function(x){as.matrix(iris[iris$Species == x, 1:4])}
#'                  )
#'                  
#' names(iris_ls) <- irisSpecies
#' 
#' BoxesM(iris_ls)
BoxesM <- function(x, ...){
  ls <- lazy_dots(...)
  matrix_ls <- x

  statistic <- BoxesMStat(matrix_ls)


  out <- list(statistic = statistic)
  out
}

