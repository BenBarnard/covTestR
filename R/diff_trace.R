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
#' @examples diff_trace(mcSamples(c(0,0,0), diag(1, 3), 10, 2, matrix = FALSE, tidy = TRUE), group = population, variables = variable, samples = samples, value = value, tidy = TRUE)
#'
diff_trace <- function(data, ...) {
  UseMethod("diff_trace")
}

#' @export
#'
#' @importFrom lazyeval expr_find
#'
diff_trace.data.frame <- function(x, group, ..., variables, samples, value, tidy = FALSE){
  if(tidy == TRUE){
    tidyDataDftoMatrix(data = x,
                       group = expr_find(group),
                       variables = expr_find(variable),
                       samples = expr_find(samples),
                       value = expr_find(value),
                       test = expr_find(diff_trace.matrix))
  }else{
    dataDftoMatrix(data = x,
                   group = expr_find(group),
                   test = expr_find(diff_trace.matrix))
  }
}

#' @export
#'
#' @importFrom plyr llply
#' @importFrom plyr mlply
#'
diff_trace.matrix<- function(...){
  matrix_ls <- list(...)

  A_ls <- llply(matrix_ls, A_func)

  n <- llply(matrix_ls, function(matrix){
    nrow(matrix)
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

  trace1 - trace2
}
