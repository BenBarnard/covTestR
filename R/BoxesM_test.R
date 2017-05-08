#' Test of Equality of Covariances given by Box's M
#'
#' @inheritParams Chaipitak2013
#'
#' @return Test statistic for Chaipitak 2013
#'
#' @export
#'
#' @examples BoxsM_test(iris, group = Species)
#'
BoxesM_test <- function(x, ...){
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

  n_overall <- Reduce(`+`, lapply(ns, function(x){x - 1}))

  statistic <- n_overall * log(det(overall_cov)) - Reduce(`+`, mapply(function(x, y){x * log(det(y))}, ns, sample_covs))
}

