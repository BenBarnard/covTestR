




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
#' @return Test statistic for Srivastava 2007
#'
#' @export
#'
#' @examples Srivastava2007_test(mcSamples(c(0,0,0), diag(1, 3), 10, 2))
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

#' Test of Equality of Covariances given by Srivastava and Yanagihara 2010
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
#' @return Test statistic for Srivastava and Yanagihara 2010
#'
#' @export
#'
#' @examples SrivastavaYanagihara2010_test(mcSamples(c(0,0,0), diag(1, 3), 10, 2))
SrivastavaYanagihara2010_test <- function(data, ...){
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
    ahat2_ls <- cbind(n, p, sample_covs) %>% mlply(ahat2i_func)
    ahat1_ls <- cbind(p, sample_covs) %>% mlply(function(p, sample_covs){
      ahat1_func(p, sample_covs)
    })
    ahat3 <- ahat3_func(A_ls[[1]], A_ls[[2]], p[[1]], n[[1]], n[[2]], ahat2, ahat1)
    ahat4 <- ahat4_func(A_ls[[1]], A_ls[[2]], p[[1]], n[[1]], n[[2]], ahat2, ahat1)
    ksihat2_ls <- cbind(n, p, ahat1, ahat2, ahat3, ahat4) %>% mlply(ksihat2i_func)
    gammahat_ls <- cbind(ahat2_ls, ahat1_ls) %>% mlply(function(ahat2_ls, ahat1_ls){
      gammahati_func(ahat2_ls, ahat1_ls)
      })
  }
  gammahat1 <- gammahat_ls[[1]]
  gammahat2 <- gammahat_ls[[2]]
  ksihat21 <- ksihat2_ls[[1]]
  ksihat22 <- ksihat2_ls[[2]]
  ((gammahat1 - gammahat2) / sqrt(ksihat21 + ksihat22))
}
