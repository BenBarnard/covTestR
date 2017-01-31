#' Test of Equality of Covariances given by Srivastava et al. 2014
#'
#' Performs 2 and k sample equality of covariance matrix test using Srivastava et al. 2014
#'
#' @param x data as data.frame, grouped_df, resample or matrix object
#' @param group group or population variable
#' @param ... other options passed to functions
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
                   method = expr_find(Srivastava2014_test.matrix),
                   .dots = lazy_dots(...))
}

#' @export
#' @rdname Srivastava2014_test
#' @importFrom lazyeval expr_find
#' @importFrom lazyeval lazy_dots
Srivastava2014_test.grouped_df <- function(x, ...){
  dataDftoMatrix(data = x,
                 group = attributes(x)$vars[[1]],
                 test = expr_find(Srivastava2014_test.matrix))
}

#' @export
#' @rdname Srivastava2014_test
#' @importFrom lazyeval expr_find
#' @importFrom lazyeval lazy_dots
Srivastava2014_test.resample <- function(x, ...){
  x <- as.data.frame(x)
  dataDftoMatrix(data = x,
                 group = attributes(x)$vars[[1]],
                 method = expr_find(Srivastava2014_test.matrix),
                 .dots = lazy_dots(...))
}

#' @export
#' @rdname Srivastava2014_test
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#' @importFrom stats cov
#'
Srivastava2014_test.matrix<- function(...){
  ls <- lazy_dots(...)
  matrix_ls <- lazy_eval(ls[str_detect(names(ls), "x.")])
  names(matrix_ls) <- str_replace(names(matrix_ls), "x.", "")

  ns <- lapply(matrix_ls, function(matrix){
    nrow(matrix)
  })

  p <- lapply(matrix_ls, function(matrix){
    ncol(matrix)
  })

  A_ls <- lapply(matrix_ls, A_func)

  sample_covs <- lapply(matrix_ls, cov)

  overall_cov <- overall_cov_func(A_ls, ns)

  D_ls <- lapply(matrix_ls, Di_func)

  ahat2i <- mapply(ahat2iSrivastava2014_func, ns, p, D_ls, A_ls, SIMPLIFY = FALSE)

  ahat2 <- ahat2Srivastava2014_func(ahat2i, ns)

  Srivastava2014(ns, p[[1]], ahat2, ahat2i, sample_covs)
}

#' @keywords internal
Srivastava2014 <- function(ns, p, ahat2, ahat2i, sample_covs){
  comb <- combn(length(ns), 2, simplify = FALSE)
  theta <- 4 * (ahat2 ^ 2) * (Reduce(`+`, lapply(comb, function(x){
    ((1 / ns[[x[1]]]) + (1 / ns[[x[2]]])) ^ 2})) +
      (length(ns) - 1) * (length(ns) - 2) * Reduce(`+`, lapply(ns, function(x){x ^ 2})))

  Reduce(`+`, lapply(comb, function(x){
    ((ahat2i[[x[1]]] + ahat2i[[x[2]]] -
        (2 / p) * tr(sample_covs[[x[1]]] %*% sample_covs[[x[2]]])) ^ 2) / theta
  }))
}
