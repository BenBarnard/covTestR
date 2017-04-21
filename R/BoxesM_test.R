source("R/helper_functions.R")
#' Test of Equality of Covariances given by Box's M
#'
#' @inheritParams Chaipitak2013_test
#'
#' @return Test statistic for Chaipitak 2013
#'
#' @export
#'
#' @examples BoxsM_test(iris, group = Species)
#'
BoxsM_test <- function(x, ...){
  UseMethod("BoxsM_test")
}

#' @export
#' @keywords internal
BoxsM_test.list <- function(x, ...){

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

  BoxsM(ns, n_overall, sample_covs, overall_cov)
}


#' Hidden Test
#' @keywords internal
#' @param ns ns
#' @param n_overall n overall
#' @param sample_covs sample covs
#' @param overall_cov overall cov
#'
BoxsM <- function(ns, n_overall, sample_covs, overall_cov){
  n_overall * log(det(overall_cov)) - Reduce(`+`, mapply(function(x, y){x * log(det(y))}, ns, sample_covs))
}

#' @export
#' @keywords internal
BoxsM_test.data.frame <- helper(BoxsM_test)

#' @export
#' @keywords internal
BoxsM_test.grouped_df <- helper(BoxsM_test)

#' @export
#' @keywords internal
BoxsM_test.resample <- helper(BoxsM_test)

