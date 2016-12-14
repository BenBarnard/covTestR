#' Test of Equality of Covariances given by Chaipitak 2013
#'
#' @param data tidy data frame
#' @param ... other
#'
#' @return Test statistic for Chaipitak 2013
#'
#' @export
#'
#' @examples Chaipitak2013_test(mcSamples(rep(0, 100), diag(1, 100), 10, 3), group = population)
#'
Chaipitak2013_test <- function(data, ...){
  UseMethod("Chaipitak2013_test")
}

#' @export
#'
#' @importFrom lazyeval expr_find
#'
Chaipitak2013_test.data.frame <- function(x, group, ...){
  dataDftoMatrix(data = x,
                 group = expr_find(group),
                 test = expr_find(Chaipitak2013_test.matrix))
}

#' @export
#'
#' @importFrom lazyeval expr_find
#'
Chaipitak2013_test.grouped_df <- function(x, ...){
  dataDftoMatrix(data = x,
                 group = attributes(x)$vars[[1]],
                 test = expr_find(Chaipitak2013_test.matrix))
}

#' @export
#'
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_replace
#' @importFrom stringr str_detect
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

  Chaipitak2013_test.default(ahat2i, deltahat2)
}

#' @export
#'
Chaipitak2013_test.default <- function(ahat2i, deltahat2){
  comb <- combn(length(ahat2i), 2, simplify = FALSE)
  Reduce(`+`, lapply(comb, function(x){
    (max(sapply(ahat2i[x], function(x){x})) / min(sapply(ahat2i[x], function(x){x})) - 1) ^ 2 / deltahat2
  }))
}
