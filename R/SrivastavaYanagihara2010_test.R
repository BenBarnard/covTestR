#' Test of Equality of Covariances given by Srivastava and Yanagihara 2010
#'
#' @param data tidy data frame
#' @param ... other
#'
#' @return Test statistic for Srivastava and Yanagihara 2010
#'
#' @export
#'
#' @examples SrivastavaYanagihara2010_test(mcSamples(c(0,0,0), diag(1, 3), 10, 2), group = population)
#'
SrivastavaYanagihara2010_test <- function(data, ...){
  UseMethod("SrivastavaYanagihara2010_test")
}

#' @export
#'
#' @importFrom lazyeval expr_find
#'
SrivastavaYanagihara2010_test.data.frame <- function(x, group, ...){
    dataDftoMatrix(data = x,
                   group = expr_find(group),
                   test = expr_find(SrivastavaYanagihara2010_test.matrix))
}

#' @export
#'
#' @importFrom plyr llply
#' @importFrom plyr mlply
#'
SrivastavaYanagihara2010_test.matrix <- function(...){
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

  ahat2_ls <- mlply(cbind(n, p, sample_covs), ahat2i_func)

  ahat1_ls <- mlply(cbind(p, sample_covs), function(p, sample_covs){
    ahat1_func(p, sample_covs)
  })

  ahat3 <- ahat3_func(A_ls[[1]], A_ls[[2]], p[[1]], n[[1]], n[[2]], ahat2, ahat1)
  ahat4 <- ahat4_func(A_ls[[1]], A_ls[[2]], p[[1]], n[[1]], n[[2]], ahat2, ahat1)

  ksihat2_ls <- mlply(cbind(n, p, ahat1, ahat2, ahat3, ahat4), ksihat2i_func)

  gammahat_ls <- mlply(cbind(ahat2_ls, ahat1_ls), function(ahat2_ls, ahat1_ls){
    gammahati_func(ahat2_ls, ahat1_ls)
  })

  SrivastavaYanagihara2010_test.default(gammahat_ls[[1]], gammahat_ls[[2]], ksihat2_ls[[1]], ksihat2_ls[[2]])
}

#' @export
#'
SrivastavaYanagihara2010_test.default <- function(gammahat1, gammahat2, ksihat21, ksihat22){
  ((gammahat1 - gammahat2) / sqrt(ksihat21 + ksihat22))
}
