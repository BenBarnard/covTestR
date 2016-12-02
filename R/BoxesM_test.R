#' Test of Equality of Covariances given by Chaipitak 2013
#'
#' @param data tidy data frame
#' @param ... other
#'
#' @return Test statistic for Chaipitak 2013
#'
#' @export
#'
#' @examples BoxesM_test(mcSamples(rep(0, 10), diag(1, 10), 100, 3), group = population)
#'
BoxesM_test <- function(data, ...){
  UseMethod("BoxesM_test")
}

#' @export
#'
#' @importFrom lazyeval expr_find
#'
BoxesM_test.data.frame <- function(x, group, ...){
  dataDftoMatrix(data = x,
                 group = expr_find(group),
                 test = expr_find(BoxesM_test.matrix))
}

#' @export
#'
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_replace
#' @importFrom stringr str_detect
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

  BoxesM_test.default(ns, n_overall, sample_covs, overall_cov)
}

#' @export
#'
BoxesM_test.default <- function(ns, n_overall, sample_covs, overall_cov){
  n_overall * log(det(overall_cov)) - Reduce(`+`, mapply(function(x, y){x * log(det(y))}, ns, sample_covs))
}
