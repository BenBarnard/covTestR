#' Test of Equality of Covariances given by Srivastava 2007
#'
#' @param data tidy data frame
#' @param ... other
#'
#' @return Test statistic for Srivastava 2007
#'
#' @export
#'
#' @examples Srivastava2007_test(mcSamples(c(0,0,0), diag(1, 3), 10, 2))
#'
Srivastava2007_test <- function(data, ...){
  UseMethod("Srivastava2007_test")
}

#' @export
#'
#' @importFrom lazyeval expr_find
#'
Srivastava2007_test.data.frame <- function(x, group, ..., variables, samples, value, tidy = FALSE){
  if(tidy == TRUE){
    tidyDataDftoMatrix(data = x,
                       group = expr_find(group),
                       variables = expr_find(variable),
                       samples = expr_find(samples),
                       value = expr_find(value),
                       test = expr_find(Srivastava2007_test.matrix))
  }else{
    dataDftoMatrix(data = x,
                   group = expr_find(group),
                   test = expr_find(Srivastava2007_test.matrix))
  }
}

#' @export
#'
#' @importFrom plyr llply
#' @importFrom plyr mlply
#'
Srivastava2007_test.matrix <- function(...){
  matrix_ls <- list(...)

  n <- llply(matrix_ls, function(matrix){
    nrow(matrix)
  })

  p <- llply(matrix_ls, function(matrix){
    ncol(matrix)
  })

  A_ls <- llply(matrix_ls, A_func)

  sample_covs <- mlply(cbind(A_ls, n), function(A_ls, n){
    A_ls / (n - 1)
  })

  overall_cov <- (1 / (n[[1]] + n[[2]] - 2)) * (A_ls[[1]] + A_ls[[2]])
  ahat2 <- ahat2_func(n[[1]], n[[2]], p[[1]], overall_cov)
  ahat1 <- ahat1_func(p[[1]], overall_cov)
  ahat2i <- mlply(cbind(n, p, sample_covs), ahat2i_func)
  ahat4 <- ahat4_func(A_ls[[1]], A_ls[[2]], p[[1]], n[[1]], n[[2]], ahat2, ahat1)
  etahat2i <- mlply(cbind(n, p, ahat4, ahat2), etahat2i_func)

  Srivastava2007_test.default(ahat2i, etahat2i)
}

#' @export
#'
Srivastava2007_test.default <- function(ahat2i, etahat2i){
  ahat21 <- ahat2i[[1]]
  ahat22 <- ahat2i[[2]]
  etahat21 <- etahat2i[[1]]
  etahat22 <- etahat2i[[2]]

  ((ahat21 - ahat22) / (sqrt(etahat21 + etahat22)))
}
