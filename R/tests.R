#' Test of Equality of Covariances given by Schott 2007
#'
#' @param data tidy data frame
#' @param ... other
#'
#' @importFrom plyr dlply
#' @importFrom plyr llply
#' @importFrom plyr mlply
#' @importFrom plyr .
#' @importFrom magrittr %>%
#'
#' @return
#' @export
#'
#' @examples
Schott2007_test <- function(data, ...){
  if(!(missing(data))){
    matrix_ls <- dlply(data, .(Group), dataDftoMatrix)
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
  }
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


#' Test of Equality of Covariances given by Chaipitak 2013
#'
#' @param data tidy data frame
#' @param ... other
#'
#' @return
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr summarise
#' @importFrom dplyr select
#' @importFrom plyr dlply
#' @importFram plyr llply
#' @importFrom plyr mlply
#' @importFrom plyr .
#'
#'
#' @examples
Chaipitak2013_test <- function(data, ...){
  if(!(missing(data))){
    matrix_ls <- dlply(data, .(Group), dataDftoMatrix)
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
    tau <- tau_func(n[[1]], n[[2]])
    ahatStar4 <- ahatStar4_func(tau, p[[1]], overall_cov, n[[1]], n[[2]])
    deltahat2 <- deltahat2_func(ahatStar4, p[[1]], ahat2, n[[1]], n[[2]])
    bhat <- bhat_func(ahat2i[[1]], ahat2i[[2]])
  }
  bhat <- bhat
  deltahat2 <- deltahat2
  (bhat - 1) / sqrt(deltahat2)
}

#' Test of Equality of Covariances given by Srivastava 2007
#'
#' @param data tidy data frame
#' @param ... other
#'
#' @importFrom plyr dlply
#' @importFrom plyr llply
#' @importFrom plyr mlply
#' @importFrom plyr .
#' @importFrom magrittr %>%
#'
#' @return
#' @export
#'
#' @examples
Srivastava2007_test <- function(data, ...){
  if(!(missing(data))){
    matrix_ls <- dlply(data, .(Group), dataDftoMatrix)
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
  }
  ahat21 <- ahat2i[[1]]
  ahat22 <- ahat2i[[2]]
  etahat21 <- etahat2i[[1]]
  etahat22 <- etahat2i[[2]]
  ((ahat21 - ahat22) / (sqrt(etahat21 + etahat22)))
}
