SrivastavaYanagihara2010 <- if("covariance" %in% class(x[[1]])){
  ns <- lapply(matrix_ls, function(matrix){
    attributes(matrix)$df + 1
  })

  p <- lapply(matrix_ls, function(matrix){
    ncol(matrix)
  })

  sample_covs <- matrix_ls

  A_ls <- mapply(function(sample_covs, ns){
    sample_covs * (ns - 1)
  }, sample_covs = sample_covs, ns = ns, SIMPLIFY = FALSE)
}
