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
#' @examples Chaipitak2013_test(iris, group = Species)
#'
Chaipitak2013_test <- function(x, ...){
  UseMethod("Chaipitak2013_test")
}

#' @export
#' @keywords internal
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_replace
#' @importFrom stringr str_detect
#' @importFrom stats cov
#' @importFrom stats pchisq
#'
Chaipitak2013_test.list <- function(x, ...){

  ls <- lazy_dots(...)
  matrix_ls <- x

  if(!("covariance" %in% class(x[[1]])) & ("matrix" %in% class(x[[1]]))){
  ns <- lapply(matrix_ls, function(matrix){
    nrow(matrix)
  })

  p <- lapply(matrix_ls, function(matrix){
    ncol(matrix)
  })

  A_ls <- lapply(matrix_ls, A_func)

  sample_covs <- lapply(matrix_ls, cov)
  }

  if("covariance" %in% class(x[[1]])){
    ns <- lapply(matrix_ls, function(matrix){
      attributes(matrix)$n
    })

    p <- lapply(matrix_ls, function(matrix){
      ncol(matrix)
    })

    sample_covs <- matrix_ls

    A_ls <- mapply(function(sample_covs, ns){
      sample_covs * (ns - 1)
    }, sample_covs = sample_covs, ns = ns, SIMPLIFY = FALSE)
  }

  overall_cov <- overall_cov_func(A_ls, ns)

  ahat2i <- mapply(ahat2i_func, ns, p, sample_covs, SIMPLIFY = FALSE)
  ahat2 <- ahat2_func(ns, overall_cov, p[[1]])

    tau <- tau_func(ns)
    ahatStar4 <- ahatStar4_func(tau, p[[1]], overall_cov, ns)
    deltahat2 <- mapply(deltahat2_func, ahatStar4, p, ahat2, ns, SIMPLIFY = FALSE)
    b <- mapply(function(ahat2, ahat2i){ahat2i / ahat2}, ahat2, ahat2i, SIMPLIFY = FALSE)

  xmin <- names(matrix_ls[1])
  xmax <- names(matrix_ls[length(matrix_ls)])
  xother <- names(matrix_ls[-c(1, length(matrix_ls))])

  data.name <- Reduce(paste0, past(xmin = xmin, xother, xmax = xmax))

  statistic <- Chaipitak2013(b, deltahat2)
  names(statistic) <- "Chi Squared"

  parameter <- length(matrix_ls) - 1
  names(parameter) <- "df"

  null.value <- 0
  names(null.value) <- "difference in covariances"

  p.value <- 1 - pchisq(statistic, parameter)

  estimate <- sample_covs
  names(estimate) <- paste0("covariance of ", names(matrix_ls))

  estimate <- if(nrow(estimate[[1]]) > 5){
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
              method = "Chaipitak and Chongchareon 2013 Equality of Covariance Test",
              data.name = data.name)
  class(obj) <- "htest"
  obj
}

#' @keywords internal
Chaipitak2013 <- function(b, deltahat2){
  Reduce(`+`, mapply(function(b, deltahat2){
    (b - 1) ^ 2 / deltahat2
  }, b, deltahat2, SIMPLIFY = FALSE))
}

#' @export
#' @keywords internal
Chaipitak2013_test.data.frame <- Chaipitak2013_test.resample <- Chaipitak2013_test.grouped_df <- helper(Chaipitak2013_test)

