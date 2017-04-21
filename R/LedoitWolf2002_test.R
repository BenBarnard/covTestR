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
LedoitWolf2002_test <- function(x, ...){
  UseMethod("LedoitWolf2002_test")
}

#' @export
#' @keywords internal
#' @importFrom stats cov
#' @importFrom stats pchisq
#'
LedoitWolf2002_test.covariance <- function(x, covMat = "Identity", ...){
  p <- ncol(x)
  n <- attributes(x)$df + 1
  S <- x
  if(covMat == "Identity"){covMat <- diag(1, p)}

  statistic <- LedoitWolf2002(n, p, S, covMat)
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
LedoitWolf2002_test.matrix <- function(x, ...){
  p <- ncol(x)
  n <- nrow(x)
  S <- cov(x)
if(!(exists("covMat"))){covMat <- diag(1, p)}

  statistic <- LedoitWolf2002(n, p, S, covMat)
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
LedoitWolf2002 <- function(n, p, S, covMat){
  inv <- solve(covMat)
  mid <- (S - covMat) %*% inv
  n * p * (tr(mid %*% mid) / p -
    (p / n) * ((tr(S) / p) ^ 2) +
    (p / n)) / 2
}

#' @export
#' @keywords internal
LedoitWolf2002_test.data.frame <- helperOne(LedoitWolf2002_test)

#' @export
#' @keywords internal
LedoitWolf2002_test.resample <- helperOne(LedoitWolf2002_test)
