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
#' @examples Ahmad2015(iris[1:50, 1:3])
#'
Ahmad2015 <- function(x, Sigma = "identity", ...){
  UseMethod("Ahmad2015")
}

#' @export
#' @keywords internal
#' @importFrom stats cov
#' @importFrom stats pchisq
#'
Ahmad2015.covariance <- function(x, Sigma = "identity", ...){
  p <- ncol(x)
  n <- attributes(x)$df + 1
  S <- x

  if(Sigma[[1]] == "identity"){
    svCov <- svd(x)
    x_ <- svCov$u %*% diag(sqrt(svCov$d))
  }else{
    svCov <- svd(x)
    sv <- svd(Sigma)
    x_ <- svCov$u %*% diag(sqrt(svCov$d)) %*%
      solve(sv$u %*% diag(sqrt(sv$d)))
  }


  statistic <- Ahmad2015_(n, p, x_)
  names(statistic) <- "Normal"

  parameter <- c(0, 4 * (2 / (p / n + 1)))
  names(parameter) <- c("Mean", "Variance")

  null.value <- 0
  names(null.value) <- "difference between the Sample Covariance and the Null Covarince Structure"

  p.value <- 1 - pnorm(abs(statistic), 0, 4 * (2 / (p / n + 1)))

  estimate <- S
  estimate <- if(nrow(estimate) > 5){
    NULL
  }else{
    estimate
  }

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
Ahmad2015.matrix <- function(x, Sigma = "identity", ...){
  p <- ncol(x)
  n <- nrow(x)
  S <- cov(x)

  if(Sigma[[1]] == "identity"){
    x_ <- x
  }else{
    sv <- svd(Sigma)
    x_ <- x %*% solve(sv$u %*% diag(sqrt(sv$d)))
  }

  statistic <- Ahmad2015_(n, p, x_)
  names(statistic) <- "Normal"

  parameter <- c(0, 4 * (2 / (p / n + 1)))
  names(parameter) <- c("Mean", "Variance")

  null.value <- 0
  names(null.value) <- "difference between the Sample Covariance and the Null Covarince Structure"

  p.value <- 1 - pnorm(abs(statistic), 0, 4 * (2 / (p / n + 1)))

  estimate <- S
  estimate <- if(nrow(estimate) > 5){
    NULL
  }else{
    estimate
  }

  obj <- list(statistic = statistic,
              parameter = parameter,
              p.value = p.value,
              estimate = estimate,
              null.value = null.value,
              alternative = "two.sided",
              method = "Nagao 1973 Test of Covariance Structure")
  class(obj) <- "htest"
  obj
}

#' @keywords internal
Ahmad2015_ <- function(n, p, x_){
  n * (c3(x_) / p - 2 * c1(x_) / p + 1)
}
