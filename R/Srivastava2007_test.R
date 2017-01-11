#' Test of Equality of Covariances given by Srivastava 2007
#'
#' @param data tidy data frame
#' @param ... other
#'
#' @return Test statistic for Srivastava 2007
#'
#' @export
#'
#' @examples Srivastava2007_test(iris, group = Species)
#'
Srivastava2007_test <- function(data, ...){
  UseMethod("Srivastava2007_test")
}

#' @export
#' @rdname Srivastava2007_test
#' @importFrom lazyeval expr_find
#'
Srivastava2007_test.data.frame <- function(x, group, ...){
  dataDftoMatrix(data = x,
                 group = expr_find(group),
                 test = expr_find(Srivastava2007_test.matrix))
}

#' @export
#' @rdname Srivastava2007_test
#' @importFrom lazyeval expr_find
#'
Srivastava2007_test.grouped_df <- function(x, ...){
  dataDftoMatrix(data = x,
                 group = attributes(x)$vars[[1]],
                 test = expr_find(Srivastava2007_test.matrix))
}

#' @export
#' @rdname Srivastava2007_test
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#'
Srivastava2007_test.matrix <- function(...){
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

  ahat1 <- ahat1_func(p[[1]], overall_cov)

  ahat4 <- ahat4_func(A_ls, p[[1]], ns, ahat2, ahat1)

  etahat2i <- mapply(etahat2i_func, ns, p, ahat4, ahat2, SIMPLIFY = FALSE)


  Srivastava2007_test.default(ahat2i, etahat2i)
}

#' @export
#' @rdname Srivastava2007_test
Srivastava2007_test.default <- function(ahat2i, etahat2i){
  ahatbar <- Reduce(`+`, mapply(function(ahat2i, etahat2i){ahat2i / etahat2i},
                                ahat2i, etahat2i, SIMPLIFY = FALSE)) /
    Reduce(`+`, lapply(etahat2i, function(x){1 / x}))
  Reduce(`+`, mapply(function(ahat2i, etahat2i){
    ((ahat2i - ahatbar) ^ 2) / etahat2i
  }, ahat2i, etahat2i, SIMPLIFY = FALSE))
}
