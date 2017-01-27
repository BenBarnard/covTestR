#' Test of Equality of Covariances given by Chaipitak and Chongcharoen 2013
#'
#' Performs 2 and k sample equality of covariance matrix test using Chaipitak and Chongcharoen 2013
#'
#' @param x data as data.frame, grouped_df, resample or matrix object
#' @param group group or population variable
#' @param ... other options passed to functions
#'
#' @name Chaipitak2013 Test
#' @section Usage:
#' \preformatted{Chaipitak2013_test(x, ...)
#'
#' ## S3 method for class 'data.frame'
#' Chaipitak2013_test(x, group, ...)
#'
#' ## S3 method for class 'grouped_df'
#' Chaipitak2013_test(x, ...)
#'
#' ## S3 method for class 'resample'
#' Chaipitak2013_test(x, ...)
#'
#' ## S3 method for class 'matrix'
#' Chaipitak2013_test(...)}
#'
#' @return Test statistic of the hypothesis test
#'
#' @details \deqn{T_{Sc} = \sum \limits^{k}_{i < j} \frac{ \left( \hat{a}_{2i_{Sc}} + \hat{a}_{2j_{Sc}} - \frac{2}{p}tr \left( S_i S_j \right) \right) ^ 2}{\theta_{Sc}}}
#'
#' @references Chaipitak, S. and Chongcharoen, S. (2013). A test for testing the equality
#' of two covariance matrices for high-dimensional data. Journal of Applied Sciences, 13(2):270–277.
#'
#' @examples Chaipitak2013_test(iris, group = Species)
NULL

#' Test of Equality of Covariances given by Chaipitak and Chongcharoen 2013
#'
#' Performs 2 and k sample equality of covariance matrix test using Chaipitak and Chongcharoen 2013
#'
#' @param x data as data.frame, grouped_df, resample or matrix object
#' @param group group or population variable
#' @param ... other options passed to functions
#'
#' @return Test statistic of the hypothesis test
#'
#' @details \deqn{T_{Sc} = \sum \limits^{k}_{i < j} \frac{ \left( \hat{a}_{2i_{Sc}} + \hat{a}_{2j_{Sc}} - \frac{2}{p}tr \left( S_i S_j \right) \right) ^ 2}{\theta_{Sc}}}
#'
#' @export
#' @keywords internal
#'
#' @references Chaipitak, S. and Chongcharoen, S. (2013). A test for testing the equality
#' of two covariance matrices for high-dimensional data. Journal of Applied Sciences, 13(2):270–277.
#'
#' @examples Chaipitak2013_test(iris, group = Species)
#'
Chaipitak2013_test <- function(x, ...){
  UseMethod("Chaipitak2013_test")
}

#' @export
#' @rdname Chaipitak2013_test
#' @importFrom lazyeval expr_find
#' @importFrom lazyeval lazy_dots
Chaipitak2013_test.data.frame <- function(x, group, ...){
  dataDftoMatrix(data = x,
                 group = expr_find(group),
                 method = expr_find(Chaipitak2013_test.matrix),
                 .dots = lazy_dots(...))
}

#' @export
#' @rdname Chaipitak2013_test
#' @importFrom lazyeval expr_find
#' @importFrom lazyeval lazy_dots
Chaipitak2013_test.grouped_df <- function(x, ...){
  dataDftoMatrix(data = x,
                 group = attributes(x)$vars[[1]],
                 method = expr_find(Chaipitak2013_test.matrix),
                 .dots = lazy_dots(...))
}

#' @export
#' @rdname Chaipitak2013_test
#' @importFrom lazyeval expr_find
#' @importFrom lazyeval lazy_dots
Chaipitak2013_test.resample <- function(x, ...){
  x <- as.data.frame(x)
  dataDftoMatrix(data = x,
                 group = attributes(x)$vars[[1]],
                 method = expr_find(Chaipitak2013_test.matrix),
                 .dots = lazy_dots(...))
}

#' @export
#' @rdname Chaipitak2013_test
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_replace
#' @importFrom stringr str_detect
#' @importFrom stats cov
#'
Chaipitak2013_test.matrix <- function(...){
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

  ahat2i <- mapply(ahat2i_func, ns, p, sample_covs, SIMPLIFY = FALSE)
  ahat2 <- ahat2_func(ns, overall_cov, p[[1]])

    tau <- tau_func(ns)
    ahatStar4 <- ahatStar4_func(tau, p[[1]], overall_cov, ns)
    deltahat2 <- deltahat2_func(ahatStar4, p[[1]], ahat2, ns)

  Chaipitak2013(ahat2i, deltahat2)
}

#' @keywords internal
Chaipitak2013 <- function(ahat2i, deltahat2){
  comb <- combn(length(ahat2i), 2, simplify = FALSE)
  Reduce(`+`, lapply(comb, function(x){
    (max(sapply(ahat2i[x], function(x){x})) / min(sapply(ahat2i[x], function(x){x})) - 1) ^ 2 / deltahat2
  }))
}
