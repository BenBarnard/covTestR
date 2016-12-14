#' Test of Equality of Covariances given by Srivastava 2014
#'
#' @param data tidy data frame
#' @param ... other
#'
#' @return Test Statistic for Srivastava 2014
#' @export
#'
#' @examples Srivastava2014_test(mcSamples(c(0,0,0), diag(1, 3), 10, 3), group = population)
#'
Srivastava2014_test <- function(data, ...) {
  UseMethod("Srivastava2014_test")
}

#' @export
#'
#' @importFrom lazyeval expr_find
#'
Srivastava2014_test.data.frame <- function(x, group, ...){
    dataDftoMatrix(data = x,
                   group = expr_find(group),
                   test = expr_find(Srivastava2014_test.matrix))
}

#' @export
#'
#' @importFrom lazyeval expr_find
#'
Srivastava2014_test.grouped_df <- function(x, ...){
  dataDftoMatrix(data = x,
                 group = attributes(x)$vars[[1]],
                 test = expr_find(Srivastava2014_test.matrix))
}

#' @export
#'
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
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

  Srivastava2014_test.default(ns, p[[1]], ahat2, ahat2i, sample_covs)
}

#' @export
#'
Srivastava2014_test.default <- function(ns, p, ahat2, ahat2i, sample_covs){
  comb <- combn(length(ns), 2, simplify = FALSE)
  theta <- 4 * (ahat2 ^ 2) * (Reduce(`+`, lapply(comb, function(x){
    ((1 / ns[[x[1]]]) + (1 / ns[[x[2]]])) ^ 2})) +
      (length(ns) - 1) * (length(ns) - 2) * Reduce(`+`, lapply(ns, function(x){x ^ 2})))

  Reduce(`+`, lapply(comb, function(x){
    ((ahat2i[[x[1]]] + ahat2i[[x[2]]] -
        (2 / p) * tr(sample_covs[[x[1]]] %*% sample_covs[[x[2]]])) ^ 2) / theta
  }))
}
