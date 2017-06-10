#' Test of Equality of Covariances given by Schott 2007
#'
#' Performs 2 and k sample equality of covariance matrix test using Schott 2007
#'
#' @inheritParams Chaipitak2013
#'
#' @return Test statistic of the hypothesis test
#'
#' @export
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#' @importFrom stats cov
#' @importFrom stats pchisq
#'
#' @references Schott, J. (2007). A test for the equality of covariance matrices when the dimension
#' is large relative to the sample sizes. Computational Statistics & Data Analysis, 51(12):6535-6542.
#'
#' @examples Schott2007(iris, group = Species)
#'
Schott2007 <- function(x, ...) {

  ls <- lazy_dots(...)
  matrix_ls <- x

  statistic <- Schott2007Stat(matrix_ls)

  xmin <- names(matrix_ls[1])
  xmax <- names(matrix_ls[length(matrix_ls)])
  xother <- names(matrix_ls[-c(1, length(matrix_ls))])

  data.name <- Reduce(paste0, past(xmin = xmin, xother, xmax = xmax))

  names(statistic) <- "Chi Squared"

  parameter <- 1
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
              method = "Schott 2007 Equality of Covariance Test",
              data.name = data.name)
  class(obj) <- "htest"
  obj
}

