#' Difference in Trace of Covariance Matrices
#'
#' @param data Data
#' @param ... other
#'
#'
#'
#' @return Difference of Traces for Covariance Matrices
#' @export
#'
#' @examples diff_trace(mcSamples(c(0,0,0), diag(1, 3), 10, 2), group = population)
#'
diff_trace <- function(data, ...) {
  UseMethod("diff_trace")
}

#' @export
#'
#' @importFrom lazyeval expr_find
#'
diff_trace.data.frame <- function(x, group, ...){
  dataDftoMatrix(data = x,
                 group = expr_find(group),
                 test = expr_find(diff_trace.matrix))
}

#' @export
#'
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_replace
#' @importFrom stringr str_detect
#'
diff_trace.matrix<- function(...){
  ls <- lazy_dots(...)
  matrix_ls <- lazy_eval(ls[str_detect(names(ls), "x.")])
  names(matrix_ls) <- str_replace(names(matrix_ls), "x.", "")

  ns <- lapply(matrix_ls, function(matrix){
    nrow(matrix)
  })

  p <- lapply(matrix_ls, function(matrix){
    ncol(matrix)
  })

  A_ls <- lapply(matrix_ls, A_func)

  sample_covs <- lapply(matrix_ls, cov)

  overall_cov <- overall_cov_func(A_ls, ns)

  ahat2i <- mapply(ahat2i_func, ns, p, sample_covs, SIMPLIFY = FALSE)

diff_trace.default(ahat2i, p[[1]], sample_covs)
}

#' @export
#'
diff_trace.default <- function(ahat2i, p, sample_covs){
  comb <- combn(length(ahat2i), 2, simplify = FALSE)

  Reduce(`+`, lapply(comb, function(x){
    (ahat2i[[x[1]]] + ahat2i[[x[2]]] -
      (2 / p) * tr(sample_covs[[x[1]]] %*% sample_covs[[x[2]]])) ^ 2
  }))
}
