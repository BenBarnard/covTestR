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
#' @importFrom plyr dlply
#' @importFrom plyr .
#' @importFrom magrittr %>%
#'
Srivastava2007_test.data.frame <- function(data, ...){
  do.call(Srivastava2007_test.matrix,
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
Srivastava2007_test.matrix <- function(...){
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
  overall_cov <- (1 / (n[[1]] + n[[2]] - 2)) * (A_ls[[1]] + A_ls[[2]])
  ahat2 <- ahat2_func(n[[1]], n[[2]], p[[1]], overall_cov)
  ahat1 <- ahat1_func(p[[1]], overall_cov)
  ahat2i <- cbind(n, p, sample_covs) %>% mlply(ahat2i_func)
  ahat4 <- ahat4_func(A_ls[[1]], A_ls[[2]], p[[1]], n[[1]], n[[2]], ahat2, ahat1)
  etahat2i <- cbind(n, p, ahat4, ahat2) %>% mlply(etahat2i_func)

  Srivastava2007_test.default(ahat2i[[1]], ahat2i[[2]], etahat2i[[1]], etahat2i[[2]])
}

#' @export
#'
Srivastava2007_test.default <- function(ahat21, ahat22, etahat21, etahat22){
  ((ahat21 - ahat22) / (sqrt(etahat21 + etahat22)))
}
