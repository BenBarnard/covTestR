#' Test of Equality of Covariances given by Srivastava 2014
#'
#' @param data tidy data frame
#' @param ... other
#'
#' @return Test Statistic for Srivastava 2014
#' @export
#'
#' @examples Srivastava2014_test(mcSamples(c(0,0,0), diag(1, 3), 10, 2))
#'
Srivastava2014_test <- function(data, ...) {
  UseMethod("Srivastava2014_test")
}

#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom plyr dlply
#' @importFrom plyr .
#'
Srivastava2014_test.data.frame <- function(data, ...){
  do.call(Srivastava2014_test.matrix,
          data %>%
            dlply(.(Group), dataDftoMatrix))
}

#' @export
#'
#' @importFrom plyr dlply
#' @importFrom plyr llply
#' @importFrom plyr mlply
#' @importFrom plyr .
#' @importFrom magrittr %>%
#'
Srivastava2014_test.matrix<- function(...){
  matrix_ls <- list(...)
  n <- matrix_ls %>% llply(function(matrix){
    nrow(matrix)
  })
  p <- matrix_ls %>% llply(function(matrix){
    ncol(matrix)
  })
  A_ls <- matrix_ls %>% llply(A_func)
  sample_covs <- cbind(A_ls, n) %>% mlply(function(A_ls, n){
    A_ls / (n - 1)
  })
  D_ls <- matrix_ls %>% llply(Di_func)
  overall_cov <- (1 / (n[[1]] + n[[2]] - 2)) * (A_ls[[1]] + A_ls[[2]])
  ahat2i <- cbind(n, p, D_ls, A_ls) %>%
    mlply(function(n, p, D_ls, A_ls){
    ahat2iSrivastava2014_func(n, p, D_ls, A_ls)
  })
  ahat2 <- ahat2Srivastava2014_func(ahat2i[[1]], ahat2i[[2]], n[[1]], n[[1]])

  Srivastava2014_test.default(n, p, ahat2, ahat2i, sample_covs)
}

#' @export
#'
Srivastava2014_test.default <- function(n, p, ahat2, ahat2i, sample_covs){
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
