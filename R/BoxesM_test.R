#' Test of Equality of Covariances given by Chaipitak 2013
#'
#' @inheritParams Chaipitak2013_test
#'
#' @return Test statistic for Chaipitak 2013
#'
#' @export
#'
#' @examples BoxesM_test(iris, group = Species)
#'
BoxesM_test <- function(x, ...){
  UseMethod("BoxesM_test")
}

#' @export
#' @rdname BoxesM_test
#' @importFrom dplyr select
BoxesM_test.data.frame <- function(x, group, ...){
  dataDftoMatrix(data = x,
                 group = lazyeval::expr_find(group),
                 method = lazyeval::expr_find(BoxesM_test.list),
                 .dots = lazyeval::lazy_dots(...))
}

#' @export
#' @rdname BoxesM_test
BoxesM_test.grouped_df <- function(x, ...){
  dataDftoMatrix(data = x,
                 group = attributes(x)$vars[[1]],
                 test = lazyeval::expr_find(BoxesM_test.list))
}

#' @export
#' @rdname BoxesM_test
BoxesM_test.list <- function(x, ...){

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

  n_overall <- Reduce(`+`, lapply(ns, function(x){x - 1}))

  BoxesM(ns, n_overall, sample_covs, overall_cov)
}


#' Hidden Test
#'
#' @param ns ns
#' @param n_overall n overall
#' @param sample_covs sample covs
#' @param overall_cov overall cov
#'
BoxesM <- function(ns, n_overall, sample_covs, overall_cov){
  n_overall * log(det(overall_cov)) - Reduce(`+`, mapply(function(x, y){x * log(det(y))}, ns, sample_covs))
}
