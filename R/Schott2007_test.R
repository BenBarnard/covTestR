#' Test of Equality of Covariances given by Schott 2007
#'
#' @param data tidy data frame
#' @param ... other
#'
#' @return Test Statistic for Schott 2007
#' @export
#'
#' @examples Schott2007_test(mcSamples(c(0,0,0), diag(1, 3), 10, 2), group = population)
#'
Schott2007_test <- function(data, ...) {
  UseMethod("Schott2007_test")
}

#' @export
#'
#' @importFrom lazyeval expr_find
#'
Schott2007_test.data.frame <- function(x, group, ...){
  dataDftoMatrix(data = x,
                 group = expr_find(group),
                 test = expr_find(Schott2007_test.matrix))
}

#' @export
#'
#' @importFrom plyr llply
#' @importFrom plyr mlply
#'
Schott2007_test.matrix<- function(...){
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

  Schott2007_test.default(n, p, ahat2, ahat2i, sample_covs)
}

#' @export
#'
Schott2007_test.default <- function(n, p, ahat2, ahat2i, sample_covs){
  n1 <- n[[1]]
  n2 <- n[[2]]
  p <- p[[1]]
  ahat2 <- ahat2
  ahat21 <- ahat2i[[1]]
  ahat22 <- ahat2i[[2]]
  sample_cov1 <- sample_covs[[1]]
  sample_cov2 <- sample_covs[[2]]

  (((n1 - 1) * (n2 - 1)) / (2 * (n1 + n2 - 2) * ahat2)) *
    (ahat21 + ahat22 - (2 / p) * tr(sample_cov1 %*% sample_cov2))
}
