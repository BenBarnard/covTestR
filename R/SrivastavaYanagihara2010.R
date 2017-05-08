source("R/helper_functions.R")
#' Test of Equality of Covariances given by Srivastava and Yanagihara 2010
#'
#' @inheritParams Chaipitak2013
#'
#' @return Test statistic for Srivastava and Yanagihara 2010
#'
#' @export
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#' @importFrom stats cov
#' @importFrom stats pchisq
#'
#' @references Srivastava, M. and Yanagihara, H. (2010). Testing the equality of several covariance matrices with
#' fewer observation that the dimension. Journal of Multivariate Analysis, 101(6):1319-1329.
#'
#' @examples SrivastavaYanagihara2010(iris, group = Species)
#'
SrivastavaYanagihara2010 <- function(x, ...){

  ls <- lazy_dots(...)
  matrix_ls <- x

  if(!("covariance" %in% class(x[[1]])) & ("matrix" %in% class(x[[1]]))){
  ns <- lapply(matrix_ls, function(matrix){
    nrow(matrix)
  })

  p <- lapply(matrix_ls, function(matrix){
    ncol(matrix)
  })

  A_ls <- lapply(matrix_ls, A_func)

  sample_covs <- lapply(matrix_ls, cov)
  }

  if("covariance" %in% class(x[[1]])){
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

  overall_cov <- overall_cov_func(A_ls, ns)

  ahat2 <- ahat2_func(ns, overall_cov, p[[1]])

  ahat1 <- ahat1_func(p[[1]], overall_cov)

  ahat2i <- mapply(ahat2i_func, ns, p, sample_covs, SIMPLIFY = FALSE)

  ahat1i <- lapply(sample_covs, function(x){ahat1i_func(p[[1]], x)})

  ahat3 <- ahat3_func(A_ls, p[[1]], ns, ahat2, ahat1)

  ahat4 <- ahat4_func(A_ls, p[[1]], ns, ahat2, ahat1)

  ksihat2_ls <- mapply(ksihat2i_func, ns, p, ahat1, ahat2,
                       ahat3, ahat4, SIMPLIFY = FALSE)

  gammahat_ls <- mapply(gammahati_func, ahat2i, ahat1i, SIMPLIFY = FALSE)

  xmin <- names(matrix_ls[1])
  xmax <- names(matrix_ls[length(matrix_ls)])
  xother <- names(matrix_ls[-c(1, length(matrix_ls))])

  data.name <- Reduce(paste0, past(xmin = xmin, xother, xmax = xmax))

  gammahatbar <- Reduce(`+`, mapply(
    function(gammahat_ls, ksihat2_ls){
      gammahat_ls / ksihat2_ls
    }, gammahat_ls, ksihat2_ls, SIMPLIFY = FALSE)) /
    Reduce(`+`, lapply(ksihat2_ls, function(x){1 / x}))

  statistic <- Reduce(`+`, mapply(function(gammahat_ls, ksihat2_ls){
    ((gammahat_ls - gammahatbar) ^ 2) / ksihat2_ls
  }, gammahat_ls, ksihat2_ls, SIMPLIFY = FALSE))

  names(statistic) <- "Chi Squared"

  parameter <- length(matrix_ls) - 1
  names(parameter) <- "df"

  null.value <- 0
  names(null.value) <- "difference in covariances"

  p.value <- 1 - pchisq(statistic, parameter)

  estimate <- sample_covs
  names(estimate) <- paste0("covariance of ", names(matrix_ls))

  estimate <- if(nrow(estimate[[1]]) > 5){
    NULL
  }else{
    estimate
  }

  obj <- list(statistic = statistic,
              parameter = parameter,
              p.value = p.value,
              estimate = estimate,
              null.value = null.value,
              alternative = "two.sided",
              method = "Srivastava and Yanagihara 2010 Equality of Covariance Test",
              data.name = data.name)
  class(obj) <- "htest"
  obj
}
