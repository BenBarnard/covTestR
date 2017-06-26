#' Test of Equality of Covariances given by Chaipitak and Chongcharoen 2013
#'
#' Performs 2 and k sample equality of covariance matrix test using Chaipitak and Chongcharoen 2013
#'
#' @inheritParams structureCovariances
#'
#' @return Test statistic of the hypothesis test
#'
#'
#' @export
#'
#' @references Chaipitak, S. and Chongcharoen, S. (2013). A test for testing the equality of two covariance
#' matrices for high-dimensional data. Journal of Applied Sciences, 13(2):270-277.
#'
#' @examples Chen2010(iris[1:50, 1:3])
#'
Chen2010 <- function(x, Sigma = "identity", ...){
  UseMethod("Chen2010")
}

#' @export
#' @keywords internal
#' @importFrom stats cov
#' @importFrom stats pnorm
#'
Chen2010.covariance <- function(x, Sigma = "identity", ...){
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


  statistic <- Chen2010_(n, p, x_)
  names(statistic) <- "Standard Normal"

  parameter <- c(0, 1)
  names(parameter) <- c("Mean", "Variance")

  null.value <- 0
  names(null.value) <- "difference between the Sample Covariance and the Null Covarince Structure"

  p.value <- 1 - pnorm(abs(statistic))

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
#' @importFrom stats pnorm
#'
Chen2010.matrix <- function(x, Sigma = "identity", ...){
  p <- ncol(x)
  n <- nrow(x)
  S <- cov(x)

  if(Sigma[[1]] == "identity"){
    x_ <- x
  }else{
    sv <- svd(Sigma)
    x_ <- x %*% solve(sv$u %*% diag(sqrt(sv$d)))
  }

  statistic <- Chen2010_(n, p, x_)
  names(statistic) <- "Standard Normal"

  parameter <- c(0, 1)
  names(parameter) <- c("Mean", "Variance")

  null.value <- 0
  names(null.value) <- "difference between the Sample Covariance and the Null Covarince Structure"

  p.value <- 1 - pnorm(abs(statistic))

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
Chen2010_ <- function(n, p, x_){
  n * (bilinearsquare(x_) / p -
         2 * bilinearcube(x_) / p +
         bilinearquad(x_) / p -
         2 * quadra(x_) / p +
         2 * bilinearoff(x_) / p +
         1) / 2
}
