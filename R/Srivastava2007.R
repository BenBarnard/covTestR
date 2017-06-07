#' Test of Equality of Covariances given by Srivastava 2007
#'
#' Performs 2 and k sample equality of covariance matrix test using Srivastava 2007
#'
#' @inheritParams Chaipitak2013
#'
#' @return Test statistic of the hypothesis test
#'
#' @export
#'
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#' @importFrom stats cov
#' @importFrom stats pchisq
#'
#' @references Srivastava, M. S. (2007). Testing the equality of two covariance matrices and
#' independence of two sub-vectors with fewer observations than the dimension. InInternational
#' Conference on Advances in InterdisciplinaryStistics and Combinatorics, University of North Carolina
#' at Greensboro, NC, USA.
#'
#' @examples Srivastava2007(iris, group = Species)
#'
Srivastava2007 <- function(x, ...){

  ls <- lazy_dots(...)
  matrix_ls <- x

  if(!("covariance" %in% class(x[[1]])) & ("matrix" %in% class(x[[1]]))){
    statistic <- Srivastava2007Stat(matrix_ls)
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


  overall_cov <- overall_cov_func(A_ls, ns)

  ahat2i <- mapply(ahat2i_func, ns, p, sample_covs, SIMPLIFY = FALSE)
  ahat2 <- ahat2_func(ns, overall_cov, p[[1]])

  ahat1 <- ahat1_func(p[[1]], overall_cov)

  ahat4 <- ahat4_func(A_ls, p[[1]], ns, ahat2, ahat1)

  etahat2i <- mapply(etahat2i_func, ns, p, ahat4, ahat2, SIMPLIFY = FALSE)

  ahatbar <- Reduce(`+`, mapply(function(ahat2i, etahat2i){ahat2i / etahat2i},
                                ahat2i, etahat2i, SIMPLIFY = FALSE)) /
    Reduce(`+`, lapply(etahat2i, function(x){1 / x}))

  statistic <- Reduce(`+`, mapply(function(ahat2i, etahat2i){
    ((ahat2i - ahatbar) ^ 2) / etahat2i
  }, ahat2i, etahat2i, SIMPLIFY = FALSE))

  }

  xmin <- names(matrix_ls[1])
  xmax <- names(matrix_ls[length(matrix_ls)])
  xother <- names(matrix_ls[-c(1, length(matrix_ls))])

  data.name <- Reduce(paste0, past(xmin = xmin, xother, xmax = xmax))

  names(statistic) <- "Chi Squared"

  parameter <- length(matrix_ls) - 1
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
              method = "Srivastava 2007 Equality of Covariance Test",
              data.name = data.name)
  class(obj) <- "htest"
  obj
}
