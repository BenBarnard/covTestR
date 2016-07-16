#' Sim Function to get data for Critical Value Data
#'
#' @param data the data
#' @param tests list of test to simulate
#'
#' @export
#'
#' @examples critical_value_sim(mcSamples(c(0,0,0), diag(1, 3), 10, 2))
critical_value_sim <- function(data, tests){
  UseMethod("critical_value_sim")
}

critical_value_sim.data.frame <- function(data, ...){
  do.call(critical_value_sim.matrix,
          data %>%
            dlply(.(Group), dataDftoMatrix))
}

critical_value_sim.matrix <- function(data, tests){
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
  ahat1 <- ahat1_func(p[[1]], overall_cov)
  ahat4 <- ahat4_func(A_ls[[1]], A_ls[[2]], p[[1]], n[[1]], n[[2]], ahat2, ahat1)
  etahat2i <- cbind(n, p, ahat4, ahat2) %>% mlply(etahat2i_func)
  tau <- tau_func(n[[1]], n[[2]])
  ahatStar4 <- ahatStar4_func(tau, p[[1]], overall_cov, n[[1]], n[[2]])
  deltahat2 <- deltahat2_func(ahatStar4, p[[1]], ahat2, n[[1]], n[[2]])
  bhat <- bhat_func(ahat2i[[1]], ahat2i[[2]])
  D_ls <- matrix_ls %>% llply(Di_func)
  overall_n <- n[[1]] + n[[2]]
  ahat2i_Srivastava2014 <- cbind(n, p, D_ls, A_ls) %>%
    mlply(function(n, p, D_ls, A_ls){
      ahat2iSrivastava2014_func(n, p, D_ls, A_ls)
    })
  ahat2_Srivastava2014 <- ahat2Srivastava2014_func(ahat2i[[1]], ahat2i[[2]], n[[1]], n[[1]])
  ahat2_ls <- ahat2i
  ahat1_ls <- cbind(p, sample_covs) %>% mlply(function(p, sample_covs){
    ahat1_func(p, sample_covs)
  })
  ahat3 <- ahat3_func(A_ls[[1]], A_ls[[2]], p[[1]], n[[1]], n[[2]], ahat2, ahat1)
  ahat4 <- ahat4_func(A_ls[[1]], A_ls[[2]], p[[1]], n[[1]], n[[2]], ahat2, ahat1)
  ksihat2_ls <- cbind(n, p, ahat1, ahat2, ahat3, ahat4) %>% mlply(ksihat2i_func)
  gammahat_ls <- cbind(ahat2_ls, ahat1_ls) %>% mlply(function(ahat2_ls, ahat1_ls){
    gammahati_func(ahat2_ls, ahat1_ls)
  })

  list(
    Schott2007 = Schott2007_test(n, p, ahat2, ahat2i, sample_covs),
    Chaipitak2013 = Chaipitak2013_test(bhat, deltahat2),
    Srivastava2007 = Srivastava2007_test(ahat2i[[1]], ahat2i[[2]],
                                         etahat2i[[1]], etahat2i[[2]]),
    Srivastava2014 = Srivastava2014_test.default(
      n, p, ahat2_Srivastava2014, ahat2i_Srivastava2014, sample_covs),
    SrivastavaYanagihara2010 = SrivastavaYanagihara2010_test.default(
      gammahat_ls[[1]], gammahat_ls[[2]], ksihat2_ls[[1]], ksihat2_ls[[2]])
  )
}
