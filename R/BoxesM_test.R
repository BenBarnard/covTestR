#' Test of Equality of Covariances given by Chaipitak 2013
#'
#' @param x tidy data frame
#' @param group group
#' @param ... other
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
#' @importFrom lazyeval expr_find
#'
BoxesM_test.data.frame <- function(x, group, ...){
  dataDftoMatrix(data = x,
                 group = expr_find(group),
                 method = expr_find(BoxesM_test.matrix),
                 .dots = lazy_dots(...))
}

#' @export
#' @rdname BoxesM_test
#' @importFrom lazyeval expr_find
#'
BoxesM_test.grouped_df <- function(x, ...){
  dataDftoMatrix(data = x,
                 group = attributes(x)$vars[[1]],
                 test = expr_find(BoxesM_test.matrix))
}

#' @export
#' @rdname BoxesM_test
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_replace
#' @importFrom stringr str_detect
#' @importFrom stats cov
#'
BoxesM_test.matrix <- function(...){
  ls <- lazy_dots(...)
  matrix_ls <- lazy_eval(ls[str_detect(names(ls), "x.")])
  names(matrix_ls) <- str_replace(names(matrix_ls), "x.", "")

  ns <- lapply(matrix_ls, function(matrix){
    nrow(matrix)
  })

  A_ls <- lapply(matrix_ls, A_func)

  n_overall <- Reduce(`+`, lapply(ns, function(x){x - 1}))

  sample_covs <- lapply(matrix_ls, cov)

  overall_cov <- overall_cov_func(A_ls, ns)

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
