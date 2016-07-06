#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples critical_value_sim(mcSamples(c(0,0,0), diag(1, 3), 10, 2))
critical_value_sim <- function(data){
  matrix_ls <- data %>%
    dlply(.(Group), dataDftoMatrix)
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
  ahat2i <- cbind(n, p, sample_covs) %>% mlply(ahat2i_func)

  Schott2007_test(n, p, ahat2, ahat2i, sample_covs)
}
