#' Test of Equality of Covariances given by Chaipitak 2013
#'
#' @param data tidy data frame
#' @param ... other
#'
#' @return Test statistic for Chaipitak 2013
#'
#' @export
#'
#' @examples Chaipitak2013_test(mcSamples(c(0,0,0), diag(1, 3), 10, 2, matrix = FALSE, tidy = TRUE), group = population, variables = variable, samples = samples, value = value, tidy = TRUE)
#'
Chaipitak2013_test <- function(data, ...){
  UseMethod("Chaipitak2013_test")
}

#' @export
#'
#' @importFrom lazyeval expr_find
#'
Chaipitak2013_test.data.frame <- function(x, group, ..., variables, samples, value, tidy = FALSE){
  if(tidy == TRUE){
    tidyDataDftoMatrix(data = x,
                       group = expr_find(group),
                       variables = expr_find(variable),
                       samples = expr_find(samples),
                       value = expr_find(value),
                       test = expr_find(Chaipitak2013_test.matrix))
  }else{
    dataDftoMatrix(data = x,
                   group = expr_find(group),
                   test = expr_find(Chaipitak2013_test.matrix))
  }
}

#' @export
#'
#' @importFrom plyr llply
#' @importFrom plyr mlply
#'
Chaipitak2013_test.matrix <- function(...){
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
    ahat2i <- mlply(cbind(n, p, sample_covs), ahat2i_func)
    tau <- tau_func(n[[1]], n[[2]])
    ahatStar4 <- ahatStar4_func(tau, p[[1]], overall_cov, n[[1]], n[[2]])
    deltahat2 <- deltahat2_func(ahatStar4, p[[1]], ahat2, n[[1]], n[[2]])
    bhat <- bhat_func(ahat2i[[1]], ahat2i[[2]])
  Chaipitak2013_test.default(bhat, deltahat2)
}

#' @export
#'
Chaipitak2013_test.default <- function(bhat, deltahat2){
  (bhat - 1) / sqrt(deltahat2)
}
