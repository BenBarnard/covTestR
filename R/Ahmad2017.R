#' Test of Homogeneity of Covariance Matrices given by Ahmad 2017
#'
#' Performs 2 and k sample homogeneity of covariance matrices test using Ahmaad 
#' 2017.
#'
#' @param x data as data.frame, grouped_df, resample or matrix object
#' @param ... other options passed to functions
#'
#' @return Test statistic of the hypothesis test
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
#' @references Ahmad, M. R. (2017). Location-invariant tests of homogeneity of 
#' large-dimensional covariance matrices. Journal of Statistical Theory and 
#' Practice, 1–15.
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
  names(null.value) <- "difference in covariances"

  p.value <- 1 - pnorm(abs(statistic))

  obj <- list(statistic = statistic,
              parameter = parameter,
              p.value = p.value,
              estimate = NULL,
              null.value = null.value,
              alternative = "two.sided",
              method = "Ahmad 2017 Homogeneity of Covariance Matrix Test",
              data.name = data.name)
  class(obj) <- "htest"
  obj
}
