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
  A_ls <- lapply(matrix_ls, A_func)

  ns <- lapply(matrix_ls, nrow)

  dualcovs <- lapply(matrix_ls, function(x){
    n <- nrow(x)
    t(t(x) - rowMeans(t(x))) %*% (t(x) - rowMeans(t(x))) / (n - 1)
  })

  dfmat <- matrix_ls
  }

  if("covariance" %in% class(x[[1]])){
    covs <- matrix_ls
    ns <- lapply(matrix_ls, function(x){attributes(x)$n})
    dfmat <- lapply(matrix_ls, function(x){
      n <- attributes(x)$n
      sv <- svd(x)
      sqdi <- diag(sqrt(sv$d))
      sv$u %*% sqdi
    })
    A_ls <- lapply(dfmat, A_func)
    dualcovs <- lapply(dfmat, function(x){
      n <- nrow(x)
      t(t(x) - rowMeans(t(x))) %*% (t(x) - rowMeans(t(x))) / (n - 1)
    })
  }
  overall_dfmat <- Reduce(rbind, dfmat)

  overall_dualcov <- t(t(overall_dfmat) - rowMeans(t(overall_dfmat))) %*%
    (t(overall_dfmat) - rowMeans(t(overall_dfmat))) /
    (nrow(overall_dfmat) - length(matrix_ls))
  overall_cov <- overall_cov_func(A_ls, ns)

  lambdahat <- lapply(covs, function(x){
    eigen(x)$values
  })

  overall_lambdahat <- eigen(overall_cov)$values

  lambdatildes <- mapply(lambdatilde, lambdahat, dualcovs, SIMPLIFY = FALSE)

  overall_lambdatilde <- lambdatilde(overall_lambdahat, overall_dualcov)

  eigdual <- lapply(dualcovs, function(x){
    eigen(x)$vectors
  })

  overall_eigdual <- eigen(overall_dualcov)$vectors



  htilde <- mapply(function(x, y, z){
    sapply(1:min(length(x), nrow(z)), function(k){

      (((nrow(z) - 1) * x[[k]]) ^ (-1 / 2)) * (t(y) - rowMeans(t(y))) %*% t(z)[,k]
    })
  }, lambdatildes, dfmat, eigdual, SIMPLIFY = FALSE)

  overall_htilde <- sapply(1:min(length(overall_lambdatilde), nrow(overall_eigdual)), function(k){

    (((nrow(overall_eigdual) - 1) * overall_lambdatilde[[k]]) ^ (-1 / 2)) * (t(overall_dfmat) - rowMeans(t(overall_dfmat))) %*% t(overall_eigdual)[,k]
  })

  ki <- mapply(function(x,y){
    tr(x) - y[[1]]
  }, dualcovs, lambdatildes, SIMPLIFY = FALSE)

  overall_ki <- tr(overall_dualcov) - overall_lambdatilde[[1]]

  xmin <- names(matrix_ls[1])
  xmax <- names(matrix_ls[length(matrix_ls)])
  xother <- names(matrix_ls[-c(1, length(matrix_ls))])

  data.name <- Reduce(paste0, past(xmin = xmin, xother, xmax = xmax))

  statistic <- Ishii2016(lambdatildes, htilde, ki, list(overall_lambdatilde), list(overall_htilde), list(overall_ki))
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
Ishii2016 <- function(lambdatildes, htilde, ki, overall_lambdatilde, overall_htilde, overall_ki){
  Reduce(`+`, mapply(function(lambdatildes, htilde, ki, overall_lambdatilde, overall_htilde, overall_ki){
    abs((lambdatildes[[1]] / overall_lambdatilde[[1]]) *
      (htilde[, 1] %*% overall_htilde[, 1] ^ -1) *
      (ki / overall_ki) - 1)
  }, lambdatildes, htilde, ki, overall_lambdatilde, overall_htilde, overall_ki, SIMPLIFY = FALSE))
}
