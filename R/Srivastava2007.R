#' Test of Equality of Covariances given by Srivastava 2007
#'
#' Performs 2 and k sample equality of covariance matrix test using Srivastava 2007
#'
#' @inheritParams Chaipitak2013
#'
#' @return Test statistic of the hypothesis test
#'
#' @export
#'
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#' @importFrom stats cov
#' @importFrom stats pchisq
#'
#' @references Srivastava, M. S. (2007). Testing the equality of two covariance matrices and
#' independence of two sub-vectors with fewer observations than the dimension. InInternational
#' Conference on Advances in InterdisciplinaryStistics and Combinatorics, University of North Carolina
#' at Greensboro, NC, USA.
#'
#' @examples Srivastava2007(iris, group = Species)
#'
Srivastava2007 <- function(x, ...){

  ls <- lazy_dots(...)
  matrix_ls <- x

  statistic <- Srivastava2007Stat(matrix_ls)

  xmin <- names(matrix_ls[1])
  xmax <- names(matrix_ls[length(matrix_ls)])
  xother <- names(matrix_ls[-c(1, length(matrix_ls))])

  data.name <- Reduce(paste0, past(xmin = xmin, xother, xmax = xmax))

  names(statistic) <- "Chi Squared"

  parameter <- length(matrix_ls) - 1
  names(parameter) <- "df"

  null.value <- 0
  names(null.value) <- "difference in covariances"

  p.value <- 1 - pchisq(statistic, parameter)


  obj <- list(statistic = statistic,
              parameter = parameter,
              p.value = p.value,
              estimate = NULL,
              null.value = null.value,
              alternative = "two.sided",
              method = "Srivastava 2007 Equality of Covariance Test",
              data.name = data.name)
  class(obj) <- "htest"
  obj
}
