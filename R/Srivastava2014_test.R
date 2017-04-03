#' Test of Equality of Covariances given by Srivastava et al. 2014
#'
#' Performs 2 and k sample equality of covariance matrix test using Srivastava et al. 2014
#'
#' @inheritParams Chaipitak2013_test
#'
#' @return Test statistic of the hypothesis test
#'
#' @export
#'
#' @details Consider the test for the equality of covariance matrices with
#' \deqn{\text{H}_0: \boldsymbol{\Sigma}_1 = \boldsymbol{\Sigma}_2 = \ldots = \boldsymbol{\Sigma}_k}
#' and
#' \deqn{\text{H}_{\text{A}}: \boldsymbol{\Sigma}_1 \neq \boldsymbol{\Sigma}_2 \neq \ldots \neq \boldsymbol{\Sigma}_k,}
#' where \eqn{\boldsymbol{\Sigma_i}} are the population covariance matrix parameters. The sample covarinance matrix estimators, \eqn{\textbf{S}_1, \textbf{S}_2, \ldots, \textbf{S}_i}, are distributed singular Wishart such that \eqn{n_i\textbf{S}_i \sim W_k\left( \boldsymbol{\Sigma}_i, n_i \right)}.
#' Then,
#' \deqn{T_{Sr2014} = \sum \limits^{k}_{i < j} \frac{ \left( \hat{a}_{2i_{Sr}} + \hat{a}_{2j_{Sr}} - \frac{2}{p}tr \left( S_i S_j \right) \right) ^ 2}{\theta_{Sr}}}
#' where,
#' \deqn{\hat{a}_{2i_{Sr}} = \frac{ \left(n_i -2 \right) \left( n_i -1 \right) tr \left( \boldsymbol{V}_i^2 \right) - n \left( n - k \right) tr \left( \boldsymbol{D}^2_i \right) + tr \left( \boldsymbol{V}_i \right)^2 }{pn_i \left( n_i -1 \right) \left( n_i -2 \right) \left( n_i -3 \right) }}
#' and
#' \deqn{\hat{a}_{2_{Sr}} = \frac{\sum \limits^k_{i = 1}(n_i -1)\hat{a}_{2i_{Sr}}}{\sum \limits^k_{i = 1}(n_i -1)}.}
#' The divisor
#' \deqn{\theta_{Sr} = 4 \hat{a}_{2_{Sr}}^2 \left( \sum \limits_{i < j}^ k \left( \frac{1}{n_i} + \frac{1}{n_j} \right) + (k - 1)(k - 2) \sum \limits_{i = 1}^k n_i^{-2} \right)}
#' is a variance term for the numerator.
#'
#' @references Srivastava, M., Yanagihara, H., and Kubokawa T. (2014). Tests for covariance
#' matrices in high dimension with less sample size. Journal of Multivariate Analysis, 130:289-309.
#'
#' @examples Srivastava2014_test(iris, group = Species)
#'
Srivastava2014_test <- function(x, ...) {
  UseMethod("Srivastava2014_test")
}

#' @export
#' @rdname Srivastava2014_test
#' @importFrom lazyeval expr_find
#' @importFrom lazyeval lazy_dots
Srivastava2014_test.data.frame <- function(x, group, ...){
    dataDftoMatrix(data = x,
                   group = expr_find(group),
                   method = expr_find(Srivastava2014_test.list),
                   .dots = lazy_dots(...))
}

#' @export
#' @rdname Srivastava2014_test
#' @importFrom lazyeval expr_find
#' @importFrom lazyeval lazy_dots
Srivastava2014_test.grouped_df <- function(x, ...){
  dataDftoMatrix(data = x,
                 group = attributes(x)$vars[[1]],
                 method = expr_find(Srivastava2014_test.list),
                 .dots = lazy_dots(...))
}

#' @export
#' @rdname Srivastava2014_test
#' @importFrom lazyeval expr_find
#' @importFrom lazyeval lazy_dots
Srivastava2014_test.resample <- function(x, ...){
  dataDftoMatrix(data = as.data.frame(x),
                 group = attributes(x)$vars[[1]],
                 method = expr_find(Srivastava2014_test.list),
                 .dots = lazy_dots(...))
}

#' @export
#' @rdname Srivastava2014_test
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#' @importFrom stats cov
#' @importFrom stats pchisq
#'
Srivastava2014_test.list <- function(x, ...){

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

    dfmat <- matrix_ls
  }

  if("covariance" %in% class(x[[1]])){
    ns <- lapply(matrix_ls, function(matrix){
      attributes(matrix)$n
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
  }

  overall_cov <- overall_cov_func(A_ls, ns)

  D_ls <- lapply(dfmat, Di_func)

  ahat2i <- mapply(ahat2iSrivastava2014_func, n = ns, p = p, D = D_ls, A = A_ls, MoreArgs = list(ns = ns), SIMPLIFY = FALSE)

  ahat2 <- ahat2Srivastava2014_func(ahat2i, ns)

  theta <- lapply(ns, function(x){2 * ahat2 / (x - 1)})

  xmin <- names(matrix_ls[1])
  xmax <- names(matrix_ls[length(matrix_ls)])
  xother <- names(matrix_ls[-c(1, length(matrix_ls))])

  data.name <- Reduce(paste0, past(xmin = xmin, xother, xmax = xmax))

  statistic <- Srivastava2014(p[[1]], ahat2, ahat2i, sample_covs, list(overall_cov), theta)
  names(statistic) <- "Chi Squared"

  parameter <- 1
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
              method = "Srivastava et al. 2014 Equality of Covariance Test",
              data.name = data.name)
  class(obj) <- "htest"
  obj
}

#' @keywords internal
Srivastava2014 <-  function(p, ahat2, ahat2i, sample_covs, overall_cov, theta){

  Reduce(`+`, mapply(function(p, ahat2i, ahat2, sample_covs, overall_cov, theta){
    (ahat2i + ahat2 - (2 / p) * tr(sample_covs %*% overall_cov) ^ 2) / (theta ^ 2)
  }, p, ahat2i, ahat2, sample_covs, overall_cov, theta, SIMPLIFY = FALSE))
}
