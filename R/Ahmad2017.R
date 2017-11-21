#' Test of Homogeneity of Covariance Matrices given by Ahmad 2017
#'
#' @inherit homogeneityCovariances
#'
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_replace
#' @importFrom stringr str_detect
#' @importFrom stats cov
#' @importFrom stats pnorm
#'
#' @export
#'
#' @references Ahmad, R. (2017). Location-invariant test of homogeneity of large-dimensional covariance matrices. Journal of Statistical Theory and Practice, 11(4):731-745. \doi{10.1080/15598608.2017.1308895}
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
#' Ahmad2017(iris_ls)
Ahmad2017 <- function(x, ...){

  ls <- lazy_dots(...)
  matrix_ls <- x
  statistic <- Ahmad2017Stat(matrix_ls)

  xmin <- names(matrix_ls[1])
  xmax <- names(matrix_ls[length(matrix_ls)])
  xother <- names(matrix_ls[-c(1, length(matrix_ls))])

  data.name <- Reduce(paste0, past(xmin = xmin, xother, xmax = xmax))

  names(statistic) <- "Standard Normal"

  parameter <- c(0, 1)
  names(parameter) <- c("Mean", "Variance")

  null.value <- 0
  names(null.value) <- "difference in covariance matrices"

  p.value <- 1 - pnorm(abs(statistic))

  obj <- list(statistic = statistic,
              parameter = parameter,
              p.value = p.value,
              estimate = NULL,
              null.value = null.value,
              alternative = "two.sided",
              method = "Ahmad 2017 Homogeneity of Covariance Matrices Test",
              data.name = data.name)
  class(obj) <- "htest"
  obj
}
