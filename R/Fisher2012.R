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
#' @examples Fisher2012(as.matrix(iris[1:50, 1:4]))
#'
Fisher2012 <- function(x, Sigma = "identity", ...){
  UseMethod("Fisher2012")
}

#' @export
#' @keywords internal
#' @importFrom stats cov
#' @importFrom stats pchisq
#'
Fisher2012.covariance <- function(x, Sigma = "identity", ...){
  p <- ncol(x)
  n <- attributes(x)$df + 1
  S <- x

  if(Sigma[[1]] == "identity"){
    S_ <- x
  }else{
    svCov <- svd(x)
    sv <- svd(Sigma)
    x_ <- svCov$u %*% diag(sqrt(svCov$d)) %*%
      solve(sv$u %*% diag(sqrt(sv$d)))
    S_ <- t(x_) %*% x_
  }

  statistic <- Fisher2012_(n - 1, p, S_)
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
#' @importFrom stats pchisq
#'
Fisher2012.matrix <- function(x, Sigma = "identity", ...){
  p <- ncol(x)
  n <- nrow(x)
  S <- cov(x)

  if(Sigma[[1]] == "identity"){
    S_ <- S
  }else{
    sv <- svd(Sigma)
    svDf <- svd(S)
    x_ <- svDf$u %*% diag(sqrt(sv$d)) %*% solve(sv$u %*% diag(sqrt(sv$d)))
    S_ <- t(x_) %*% x_
  }

  statistic <- Fisher2012_(n - 1, p, S_)
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

#' @keywords internal
Fisher2012_ <- function(n, p, S_){
  c <- p / n
  ahat2 <- ((n ^ 2) / ((n - 1) * (n + 2) * p)) *
    (tr(S_ %*% S_) - (tr(S_) ^ 2) / n)
  gamma <- ((n ^ 5) * (n ^ 2 + n + 2)) /
    ((n + 1) * (n + 2) * (n + 4) * (n + 6) * (n - 1) * (n - 2) * (n - 3))
  ahat4 <- (gamma / p) * (tr(S_ %*% S_ %*% S_ %*% S_) -
                            (4 / n) * tr(S_ %*% S_ %*% S_) * tr(S_) -
                            ((2 * (n ^ 2) + 3 * n - 6) / (n * (n ^ 2 + n + 2))) *
                            (tr(S_ %*% S_) ^ 2) +
                            ((2 * (5 * n + 6)) / (n * (n ^ 2 + n + 2))) *
                            tr(S_ %*% S_) * (tr(S_) ^ 2) -
                            ((5 * n + 6) / ((n ^ 2) * (n ^ 2 + n + 2))) *
                            (tr(S_) ^ 4))
  (n / sqrt(8 * (c ^ 2 + 12 * c + 8))) * (ahat4 - 2 * ahat2 + 1)
}
