#' Test of Homogeneity of Covariance Matrices given by Srivastava 2007
#'
#' @inherit homogeneityCovariances
#'
#' @export
#'
#' @importFrom lazyeval lazy_dots
#' @importFrom stats pchisq
#'
#' @references Srivastava, M. S. (2007). Testing the equality of two covariance matrices and independence of two sub-vectors with fewer observations than the dimension. InInternational Conference on Advances in InterdisciplinaryStistics and Combinatorics, University of North Carolina at Greensboro, NC, USA.
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
#' Srivastava2007(iris_ls)
Srivastava2007 <- function(x, ...){

  ls <- lazy_dots(...)
  matrix_ls <- x

  statistic <- Srivastava2007Stat(matrix_ls)

  xmin <- names(matrix_ls[1])
  xmax <- names(matrix_ls[length(matrix_ls)])
  xother <- names(matrix_ls[-c(1, length(matrix_ls))])

  data.name <- Reduce(paste0, past(xmin = xmin, xother, xmax = xmax))

  names(statistic) <- "Chi-Squared"

  parameter <- length(matrix_ls) - 1
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
              method = "Srivastava 2007 Homogeneity of Covariance Matrices Test",
              data.name = data.name)
  class(obj) <- "htest"
  obj
}
