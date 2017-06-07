#' Test of Equality of Covariances given by Srivastava et al. 2014
#'
#' Performs 2 and k sample equality of covariance matrix test using Srivastava et al. 2014
#'
#' @inheritParams Chaipitak2013
#'
#' @return Test statistic of the hypothesis test
#'
#' @export
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#' @importFrom stats cov
#' @importFrom stats pchisq
#'
#' @references Srivastava, M., Yanagihara, H., and Kubokawa T. (2014). Tests for covariance
#' matrices in high dimension with less sample size. Journal of Multivariate Analysis, 130:289-309.
#'
#' @examples Srivastava2014(iris, group = Species)
#'
Srivastava2014 <- function(x, ...){

  ls <- lazy_dots(...)
  matrix_ls <- x

  if(!("covariance" %in% class(x[[1]])) & ("matrix" %in% class(x[[1]]))){
    statistic <- Srivastava2014Stat(matrix_ls)
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

    dfmat <- lapply(matrix_ls, function(x){
      n <- attributes(x)$n
      sv <- svd(x)
      sqdi <- diag(sqrt(sv$d))
      sv$u %*% sqdi
    })


  overall_cov <- overall_cov_func(A_ls, ns)

  D_ls <- lapply(dfmat, Di_func)

  ahat2i <- mapply(ahat2iSrivastava2014_func, n = ns, p = p, D = D_ls, A = A_ls, MoreArgs = list(ns = ns), SIMPLIFY = FALSE)

  ahat2 <- ahat2Srivastava2014_func(ahat2i, ns)

  theta <- lapply(ns, function(x){2 * ahat2 / (x - 1)})

  statistic <- Reduce(`+`, mapply(function(p, ahat2i, ahat2, sample_covs, overall_cov, theta){
    ((ahat2i + ahat2 - (2 / p) * tr(sample_covs %*% overall_cov)) ^ 2) / (theta ^ 2)
  }, p[[1]], ahat2i, ahat2, sample_covs, list(overall_cov), theta, SIMPLIFY = FALSE))
  }

  xmin <- names(matrix_ls[1])
  xmax <- names(matrix_ls[length(matrix_ls)])
  xother <- names(matrix_ls[-c(1, length(matrix_ls))])

  data.name <- Reduce(paste0, past(xmin = xmin, xother, xmax = xmax))

  names(statistic) <- "Chi Squared"

  parameter <- 1
  names(parameter) <- "df"

  null.value <- 0
  names(null.value) <- "difference in covariances"

  p.value <- 1 - pchisq(statistic, parameter)


  obj <- list(statistic = statistic,
              parameter = parameter,
              p.value = p.value,
              estimate = NULL,
              null.value = null.value,
              alternative = "two.sided",
              method = "Srivastava et al. 2014 Equality of Covariance Test",
              data.name = data.name)
  class(obj) <- "htest"
  obj
}
