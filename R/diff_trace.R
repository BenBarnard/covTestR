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
#' @importFrom plyr llply
#' @importFrom plyr mlply
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_replace
#' @importFrom stringr str_detect
#'
diff_trace.matrix<- function(...){
  ls <- lazy_dots(...)
  matrix_ls <- lazy_eval(ls[str_detect(names(ls), "x.")])
  names(matrix_ls) <- str_replace(names(matrix_ls), "x.", "")

  A_ls <- llply(matrix_ls, A_func)

  n <- llply(matrix_ls, function(matrix){
    nrow(matrix)
  })

  p <- llply(matrix_ls, function(matrix){
    ncol(matrix)
  })

   sample_traces <- mlply(cbind(A_ls, n), function(A_ls, n){
     tr(A_ls / (n - 1))
   })

diff_trace.default(sample_traces)
}

#' @export
#'
diff_trace.default <- function(traces){
   trace1 <- traces[[1]]
   trace2 <- traces[[2]]

   (trace1 - trace2) ^ 2
}
