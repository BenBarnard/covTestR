#' Test of Equality of Covariances given by Schott 2007
#'
#' @inheritParams Chaipitak2013_test
#'
#' @references Ishii, A., Yata, K., and Aoshima, M. (2016). Asymptotic properties of the first pricipal
#' component and equality tests of covariance matrices in high-dimesion, low-sample-size context. Journal
#' of Statistical Planning and Inference, 170:186-199.
#'
#' @return Test Statistic for
#' @export
#' @keywords internal
#' @examples Ishii2016_test(iris, group = Species)
Ishii2016_test <- function(x, ...) {
  UseMethod("Ishii2016_test")
}

#' @export
#' @rdname Ishii2016_test
#' @importFrom lazyeval expr_find
Ishii2016_test.data.frame <- function(x, group, ...){
  dataDftoMatrix(data = x,
                 group = expr_find(group),
                 method = expr_find(Ishii2016_test.list),
                 .dots = lazy_dots(...))

}

#' @export
#' @rdname Ishii2016_test
#' @importFrom lazyeval expr_find
Ishii2016_test.grouped_df <- function(x, ...){
  dataDftoMatrix(data = x,
                 group = attributes(x)$vars[[1]],
                 method = expr_find(Ishii2016_test.list),
                 .dots = lazy_dots(...))
}

#' @export
#' @rdname Ishii2016_test
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#' @importFrom stats cov
Ishii2016_test.list <- function(x, ...){
  ls <- lazy_dots(...)
  matrix_ls <- x

  if(!("covariance" %in% class(x[[1]])) & ("matrix" %in% class(x[[1]]))){

  covs <- lapply(matrix_ls, cov)
  dualcovs <- lapply(matrix_ls, function(x){
    n <- nrow(x)
    A_func(t(x)) / (n - 1)
  })
  dfmat <- matrix_ls
  }

  if("covariance" %in% class(x[[1]])){
    covs <- matrix_ls
    dfmat <- lapply(matrix_ls, function(x){
      n <- attributes(x)$n
      sv <- svd(x)
      sqdi <- diag(sqrt(sv$d))
      sv$u %*% sqdi
    })
    dualcovs <- lapply(dfmat, function(x){
      n <- nrow(x)
      A_func(t(x)) / (n - 1)
    })
  }

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
  }, lambdatildes, dfmat, eigdual, SIMPLIFY = FALSE)

  ki <- mapply(function(x,y){
    tr(x) - y[[1]]
  }, dualcovs, lambdatildes, SIMPLIFY = FALSE)

  xmin <- names(matrix_ls[1])
  xmax <- names(matrix_ls[length(matrix_ls)])
  xother <- names(matrix_ls[-c(1, length(matrix_ls))])

  data.name <- Reduce(paste0, past(xmin = xmin, xother, xmax = xmax))

  statistic <- Ishii2016(lambdatildes, htilde, ki)
  names(statistic) <- "F"

  parameter <- c(length(matrix_ls), 1)
  names(parameter) <- c("df1", "df2")

  null.value <- 0
  names(null.value) <- "difference in covariances"

  p.value <- 1 - pf(statistic, df1 = parameter[1], df2 = parameter[2])

  estimate <- covs
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
              method = "Ishii 2016 Equality of Covariance Test",
              data.name = data.name)
  class(obj) <- "htest"
  obj
  }

#' @keywords internal
#' @importFrom utils combn
Ishii2016 <- function(lambdatildes, htilde, ki){
  comb <- combn(length(ki), 2, simplify = FALSE)

  Reduce(`+`, lapply(comb, function(x){
    (max(sapply(lambdatildes[x], function(x){x[[1]]})) / min(sapply(lambdatildes[x], function(x){x[[1]]}))) *
      max(abs(t(htilde[[x[1]]][,1]) %*% htilde[[x[2]]][,1]),
          abs(t(htilde[[x[1]]][,1]) %*% htilde[[x[2]]][,1]) ^ (-1)) *
      max(ki[[x[1]]] / ki[[x[2]]], ki[[x[2]]] / ki[[x[1]]])
  }))
}
