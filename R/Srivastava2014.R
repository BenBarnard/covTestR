#' Test of Equality of Covariances given by Srivastava et al. 2014
#'
#' Performs 2 and k sample equality of covariance matrix test using Srivastava et al. 2014
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
#' @references Srivastava, M., Yanagihara, H., and Kubokawa T. (2014). Tests for covariance
#' matrices in high dimension with less sample size. Journal of Multivariate Analysis, 130:289-309.
#'
#' @examples Srivastava2014(iris, group = Species)
#'
Srivastava2014 <- function(x, ...){

  ls <- lazy_dots(...)
  matrix_ls <- x
  statistic <- Srivastava2014Stat(matrix_ls)

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
              method = "Srivastava et al. 2014 Equality of Covariance Test",
              data.name = data.name)
  class(obj) <- "htest"
  obj
}
