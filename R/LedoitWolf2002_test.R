source("R/helper_functions.R")
#' Test of Equality of Covariances given by Chaipitak and Chongcharoen 2013
#'
#' Performs 2 and k sample equality of covariance matrix test using Chaipitak and Chongcharoen 2013
#'
#' @param x data as data.frame, grouped_df, resample or matrix object
#' @param ... other options passed to functions
#'
#' @return Test statistic of the hypothesis test
#'
#'
#' @export
#'
#' @references Chaipitak, S. and Chongcharoen, S. (2013). A test for testing the equality of two covariance
#' matrices for high-dimensional data. Journal of Applied Sciences, 13(2):270-277.
#'
#' @examples LedoitWolf2002_test(iris[1:50, 1:3])
#'
LedoitWolf2002 <- function(x, ...){
  UseMethod("LedoitWolf2002")
}

#' @export
#' @keywords internal
#' @importFrom stats cov
#' @importFrom stats pchisq
#'
LedoitWolf2002.covariance <- function(x, Sigma, ...){
  p <- ncol(x)
  n <- attributes(x)$df + 1
  S <- x


  statistic <- LedoitWolf2002_(n, p, S, Sigma)
  names(statistic) <- "Chi Squared"

  parameter <- p * (p + 1) / 2
  names(parameter) <- "df"

  null.value <- 0
  names(null.value) <- "difference between the Sample Covariance and the Null Covarince Structure"

  p.value <- 1 - pchisq(statistic, parameter)

  estimate <- S

  obj <- list(statistic = statistic,
              parameter = parameter,
              p.value = p.value,
              estimate = estimate,
              null.value = null.value,
              alternative = "two.sided",
              method = "Ledoit and Wolf 2002 Test of Covariance Structure")
  class(obj) <- "htest"
  obj
}

#' @export
#' @keywords internal
#' @importFrom stats cov
#' @importFrom stats pchisq
#'
LedoitWolf2002.matrix <- function(x, Sigma, ...){
  p <- ncol(x)
  n <- nrow(x)
  S <- cov(x)

  statistic <- LedoitWolf2002_(n, p, S, Sigma)
  names(statistic) <- "Chi Squared"

  parameter <- p * (p + 1) / 2
  names(parameter) <- "df"

  null.value <- 0
  names(null.value) <- "difference between the Sample Covariance and the Null Covarince Structure"

  p.value <- 1 - pchisq(statistic, parameter)

  estimate <- S

  obj <- list(statistic = statistic,
              parameter = parameter,
              p.value = p.value,
              estimate = estimate,
              null.value = null.value,
              alternative = "two.sided",
              method = "Ledoit and Wolf 2002 Test of Covariance Structure")
  class(obj) <- "htest"
  obj
}

#' @keywords internal
LedoitWolf2002_ <- function(n, p, S, Sigma){
  inv <- solve(Sigma)
  mid <- (S - Sigma) %*% inv
  n * p * (tr(mid %*% mid) / p -
    (p / n) * ((tr(S) / p) ^ 2) +
    (p / n)) / 2
}
