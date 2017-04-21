source("R/helper_functions.R")
#' Test of Equality of Covariances given by Schott 2007
#'
#' Performs 2 and k sample equality of covariance matrix test using Schott 2007
#'
#' @inheritParams Chaipitak2013_test
#'
#' @return Test statistic of the hypothesis test
#'
#' @export
#'
#' @references Schott, J. (2007). A test for the equality of covariance matrices when the dimension
#' is large relative to the sample sizes. Computational Statistics & Data Analysis, 51(12):6535-6542.
#'
#' @examples Schott2007_test(iris, group = Species)
#'
Schott2007_test <- function(x, ...) {
  UseMethod("Schott2007_test")
}

#' @export
#' @keywords internal
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#' @importFrom stats cov
#' @importFrom stats pchisq
#'
Schott2007_test.list <- function(x, ...){

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
      attributes(matrix)$df + 1
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

    theta <- lapply(ns, function(x){2 * ahat2 / (x - 1)})

    xmin <- names(matrix_ls[1])
    xmax <- names(matrix_ls[length(matrix_ls)])
    xother <- names(matrix_ls[-c(1, length(matrix_ls))])

  data.name <- Reduce(paste0, past(xmin = xmin, xother, xmax = xmax))

  statistic <- Schott2007(p[[1]], ahat2, ahat2i, sample_covs, list(overall_cov), theta)
  names(statistic) <- "Chi Squared"

  parameter <- 1
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
              method = "Schott 2007 Equality of Covariance Test",
              data.name = data.name)
  class(obj) <- "htest"
  obj
}

#' @keywords internal
#' @importFrom utils combn
Schott2007 <- function(p, ahat2, ahat2i, sample_covs, overall_cov, theta){

  Reduce(`+`, mapply(function(p, ahat2i, ahat2, sample_covs, overall_cov, theta){
    (ahat2i + ahat2 - (2 / p) * tr(sample_covs %*% overall_cov) ^ 2) / (theta ^ 2)
    }, p, ahat2i, ahat2, sample_covs, overall_cov, theta, SIMPLIFY = FALSE))
}

#' @export
#' @keywords internal
Schott2007_test.data.frame <- helper(Schott2007_test)

#' @export
#' @keywords internal
Schott2007_test.grouped_df <- helper(Schott2007_test)

#' @export
#' @keywords internal
Schott2007_test.resample <- helper(Schott2007_test)

