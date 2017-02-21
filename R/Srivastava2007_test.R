#' Test of Equality of Covariances given by Srivastava 2007
#'
#' Performs 2 and k sample equality of covariance matrix test using Srivastava 2007
#'
#' @inheritParams Chaipitak2013_test
#'
#' @return Test statistic of the hypothesis test
#'
#' @export
#'
#' @references Srivastava, M. S. (2007). Testing the equality of two covariance matrices and
#' independence of two sub-vectorswith fewer observations than the dimension. InInternational
#' Conference on Advances in InterdisciplinaryStistics and Combinatorics, University of North Carolina
#' at Greensboro, NC, USA.
#'
#' @examples Srivastava2007_test(iris, group = Species)
#'
Srivastava2007_test <- function(x, ...){
  UseMethod("Srivastava2007_test")
}

#' @export
#' @rdname Srivastava2007_test
#' @importFrom lazyeval expr_find
#'
Srivastava2007_test.data.frame <- function(x, group, ...){
  dataDftoMatrix(data = x,
                 group = expr_find(group),
                 method = expr_find(Srivastava2007_test.matrix),
                 .dots = lazy_dots(...))
}

#' @export
#' @rdname Srivastava2007_test
#' @importFrom lazyeval expr_find
#' @importFrom lazyeval lazy_dots
Srivastava2007_test.grouped_df <- function(x, ...){
  dataDftoMatrix(data = x,
                 group = attributes(x)$vars[[1]],
                 method = expr_find(Srivastava2007_test.matrix),
                 .dots = lazy_dots(...))
}

#' @export
#' @rdname Srivastava2007_test
#' @importFrom lazyeval expr_find
#' @importFrom lazyeval lazy_dots
Srivastava2007_test.resample <- function(x, ...){
  dataDftoMatrix(data = as.data.frame(x),
                 group = attributes(x)$vars[[1]],
                 method = expr_find(Srivastava2007_test.matrix),
                 .dots = lazy_dots(...))
}

#' @export
#' @rdname Srivastava2007_test
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#' @importFrom stats cov
#' @importFrom stats pchisq
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

  xmin <- names(matrix_ls[1])
  xmax <- names(matrix_ls[length(matrix_ls)])
  xother <- names(matrix_ls[-c(1, length(matrix_ls))])

  data.name <- Reduce(paste0, past(xmin = xmin, xother, xmax = xmax))

  statistic <- Srivastava2007(ahat2i, etahat2i)
  names(statistic) <- "Chi Squared"

  parameter <- length(matrix_ls) - 1
  names(parameter) <- "df"

  null.value <- 0
  names(null.value) <- "difference in covariances"

  p.value <- 1 - pchisq(statistic, parameter)

  estimate <- sample_covs
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
              method = "Srivastava 2007 Equality of Covariance Test",
              data.name = data.name)
  class(obj) <- "htest"
  obj
}



#' @keywords internal
Srivastava2007 <- function(ahat2i, etahat2i){
  ahatbar <- Reduce(`+`, mapply(function(ahat2i, etahat2i){ahat2i / etahat2i},
                                ahat2i, etahat2i, SIMPLIFY = FALSE)) /
    Reduce(`+`, lapply(etahat2i, function(x){1 / x}))
  Reduce(`+`, mapply(function(ahat2i, etahat2i){
    ((ahat2i - ahatbar) ^ 2) / etahat2i
  }, ahat2i, etahat2i, SIMPLIFY = FALSE))
}
