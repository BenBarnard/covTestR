#' Test of Equality of Covariances given by Schott 2007
#'
#' @param data tidy data frame
#' @param ... other
#'
#'
#'
#' @return Test Statistic for Schott 2007
#' @export
#'
#' @examples Ishii2016_test(mcSamples(c(0,0,0), diag(1, 3), 10, 2), group = population)
#'
Ishii2016_test <- function(data, ...) {
  UseMethod("Ishii2016_test")
}

#' @export
#'
#' @importFrom lazyeval expr_find
#'
Ishii2016_test.data.frame <- function(x, group, ...){
  dataDftoMatrix(data = x,
                 group = expr_find(group),
                 test = expr_find(Ishii2016_test.matrix))

}

#' @export
#'
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#'
Ishii2016_test.matrix<- function(...){
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

  lambdaTildes <- list(max(lambdatildes$`1`[[1]], lambdatildes$`2`[[1]]),
                       min(lambdatildes$`1`[[1]], lambdatildes$`2`[[1]]))

  hTilde <- if(max(lambdatildes$`1`[[1]], lambdatildes$`2`[[1]]) == lambdatildes$`1`[[1]]){
    max(abs(t(htilde[[1]][,1]) %*% htilde[[2]][,1]),
        abs(t(htilde[[1]][,1]) %*% htilde[[2]][,1]) ^ (-1))
    }else{
      min(abs(t(htilde[[1]][,1]) %*% htilde[[2]][,1]),
          abs(t(htilde[[1]][,1]) %*% htilde[[2]][,1]) ^ (-1))
    }

  gammaTilde <- if(max(lambdatildes$`1`[[1]], lambdatildes$`2`[[1]]) == lambdatildes$`1`[[1]]){
    max(ki[[1]] / ki[[2]], ki[[2]] / ki[[1]])
  }else{
    min(ki[[1]] / ki[[2]], ki[[2]] / ki[[1]])
  }
  Ishii2016_test.default(lambdaTildes, gammaTilde, hTilde)
}

#' @export
#'
Ishii2016_test.default <- function(lambdaTilde, gammaTilde, hTilde){
  lambdaTilde1 <- lambdaTilde[[1]]
  lambdaTilde2 <- lambdaTilde[[2]]
  browser()
  (lambdaTilde1 / lambdaTilde2) * hTilde * gammaTilde
}
