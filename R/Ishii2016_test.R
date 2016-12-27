#' Test of Equality of Covariances given by Schott 2007
#'
#' @param x data
#' @param group gropuing variable
#' @param ... other
#'
#' @return Test Statistic for
#' @export
#'
#' @examples Ishii2016_test(mcSamples(rep(0, 100), diag(1, 100), 10, 3), group = population)
Ishii2016_test <- function(data, ...) {
  UseMethod("Ishii2016_test")
}

#' @export
#' @rdname Ishii2016_test
#' @importFrom lazyeval expr_find
Ishii2016_test.data.frame <- function(x, group, ...){
  dataDftoMatrix(data = x,
                 group = expr_find(group),
                 method = expr_find(Ishii2016_test.matrix),
                 .dots = lazy_dots(...))

}

#' @export
#' @rdname Ishii2016_test
#' @importFrom lazyeval expr_find
Ishii2016_test.grouped_df <- function(x, ...){
  dataDftoMatrix(data = x,
                 group = attributes(x)$vars[[1]],
                 method = expr_find(Ishii2016_test.matrix),
                 .dots = lazy_dots(...))
}

#' @export
#' @rdname Ishii2016_test
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
Ishii2016_test.matrix <- function(...){
  ls <- lazy_dots(...)
  matrix_ls <- lazy_eval(ls[str_detect(names(ls), "x.")])
  names(matrix_ls) <- str_replace(names(matrix_ls), "x.", "")

  covs <- lapply(matrix_ls, cov)
  dualcovs <- lapply(matrix_ls, function(x){
    n <- nrow(x)
    A_func(t(x)) / (n - 1)
  })

  lambdahat <- lapply(covs, function(x){
    eigen(x)$values
  })

  lambdatildes <- mapply(lambdatilde, lambdahat, dualcovs, SIMPLIFY = FALSE)

  eigdual <- lapply(dualcovs, function(x){
    eigen(x)$vectors
  })

  htilde <- mapply(function(x, y, z){
    sapply(1:min(length(x), nrow(z)), function(k){
      (((nrow(z) - 1) * x[[k]]) ^ (-1 / 2)) * (t(y) - rowMeans(t(y))) %*% t(z)[,k]
    })
  }, lambdatildes, matrix_ls, eigdual, SIMPLIFY = FALSE)

  ki <- mapply(function(x,y){
    tr(x) - y[[1]]
  }, dualcovs, lambdatildes, SIMPLIFY = FALSE)

  Ishii2016_test.default(lambdatildes, htilde, ki)
}

#' @export
#' @rdname Ishii2016_test
Ishii2016_test.default <- function(lambdatildes, htilde, ki){
  comb <- combn(length(ki), 2, simplify = FALSE)

  Reduce(`*`, lapply(comb, function(x){
    (max(sapply(lambdatildes[x], function(x){x[[1]]})) / min(sapply(lambdatildes[x], function(x){x[[1]]}))) *
      max(abs(t(htilde[[x[1]]][,1]) %*% htilde[[x[2]]][,1]),
          abs(t(htilde[[x[1]]][,1]) %*% htilde[[x[2]]][,1]) ^ (-1)) *
      max(ki[[x[1]]] / ki[[x[2]]], ki[[x[2]]] / ki[[x[1]]])
  }))
}
