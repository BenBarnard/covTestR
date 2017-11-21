#' Test of Homogeneity of Covariance Matrices Box's M
#'
#' @inherit homogeneityCovariances
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
  
  xmin <- names(matrix_ls[1])
  xmax <- names(matrix_ls[length(matrix_ls)])
  xother <- names(matrix_ls[-c(1, length(matrix_ls))])
  
  data.name <- Reduce(paste0, past(xmin = xmin, xother, xmax = xmax))
  
  
  
  names(statistic) <- "Chi-Squared"
  
  p <- nrow(matrix_ls[[1]])
  parameter <- (length(matrix_ls) - 1) * p * (p + 1) / 2 
  names(parameter) <- "df"
  
  null.value <- 0
  names(null.value) <- "difference in covariance matrices"
  
  p.value <- 1 - pchisq(statistic, parameter)
  
  obj <- list(statistic = statistic,
              parameter = parameter,
              p.value = p.value,
              estimate = NULL,
              null.value = null.value,
              alternative = "two.sided",
              method = "Boxes' M Homogeneity of Covariance Matrices Test",
              data.name = data.name)
  class(obj) <- "htest"
  obj
}

