#' Test of Equality of Covariances given by Srivastava and Yanagihara 2010
#'
#' @param x tidy data frame
#' @param group group
#' @param ... other
#'
#' @return Test statistic for Srivastava and Yanagihara 2010
#'
#' @export
#'
#' @references Srivastava, M. and Yanagihara, H. (2010). Testing the equality of several covariance matrices with
#' fewer observation that the dimension. Journal of Multivariate Analysis, 101(6):1319-1329.
#'
#' @examples SrivastavaYanagihara2010_test(iris, group = Species)
#'
SrivastavaYanagihara2010_test <- function(x, ...){
  UseMethod("SrivastavaYanagihara2010_test")
}

#' @export
#' @rdname SrivastavaYanagihara2010_test
#' @importFrom lazyeval expr_find
#'
SrivastavaYanagihara2010_test.data.frame <- function(x, group, ...){
    dataDftoMatrix(data = x,
                   group = expr_find(group),
                   method = expr_find(SrivastavaYanagihara2010_test.matrix),
                   .dots = lazy_dots(...))
}

#' @export
#' @rdname SrivastavaYanagihara2010_test
#' @importFrom lazyeval expr_find
#'
SrivastavaYanagihara2010_test.grouped_df <- function(x, ...){
  dataDftoMatrix(data = x,
                 group = attributes(x)$vars[[1]],
                 test = expr_find(SrivastavaYanagihara2010_test.matrix))
}

#' @export
#' @rdname SrivastavaYanagihara2010_test
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#' @importFrom stats cov
#'
SrivastavaYanagihara2010_test.matrix <- function(...){
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

  ahat2 <- ahat2_func(ns, overall_cov, p[[1]])
  ahat1 <- ahat1_func(p[[1]], overall_cov)

  ahat2i <- mapply(ahat2i_func, ns, p, sample_covs, SIMPLIFY = FALSE)

  ahat1i <- lapply(sample_covs, function(x){ahat1i_func(p[[1]], x)})

  ahat3 <- ahat3_func(A_ls, p[[1]], ns, ahat2, ahat1)
  ahat4 <- ahat4_func(A_ls, p[[1]], ns, ahat2, ahat1)

  ksihat2_ls <- mapply(ksihat2i_func, ns, p, ahat1, ahat2, ahat3, ahat4, SIMPLIFY = FALSE)

  gammahat_ls <- mapply(gammahati_func, ahat2i, ahat1i, SIMPLIFY = FALSE)

  SrivastavaYanagihara2010(gammahat_ls, ksihat2_ls)
}



#' Hidden Test
#'
#' @param gammahat_ls gamma hat
#' @param ksihat2_ls ksi hat squared
#'
SrivastavaYanagihara2010 <- function(gammahat_ls, ksihat2_ls){
  gammahatbar <- Reduce(`+`, mapply(function(gammahat_ls, ksihat2_ls){gammahat_ls / ksihat2_ls},
                                    gammahat_ls, ksihat2_ls, SIMPLIFY = FALSE)) /
    Reduce(`+`, lapply(ksihat2_ls, function(x){1 / x}))
  Reduce(`+`, mapply(function(gammahat_ls, ksihat2_ls){
    ((gammahat_ls - gammahatbar) ^ 2) / ksihat2_ls
  }, gammahat_ls, ksihat2_ls, SIMPLIFY = FALSE))
}
