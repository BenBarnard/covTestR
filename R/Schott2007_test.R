#' Test of Equality of Covariances given by Schott 2007
#'
#' @param data tidy data frame
#' @param ... other
#'
#' @return Test Statistic for Schott 2007
#' @export
#'
#' @examples Schott2007_test(mcSamples(c(0,0,0), diag(1, 3), 10, 2), group = population)
#'
Schott2007_test <- function(data, ...) {
  UseMethod("Schott2007_test")
}

#' @export
#'
#' @importFrom lazyeval expr_find
#'
Schott2007_test.data.frame <- function(x, group, ...){
  dataDftoMatrix(data = x,
                 group = expr_find(group),
                 test = expr_find(Schott2007_test.matrix))
}

#' @export
#'
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#'
Schott2007_test.matrix<- function(...){
  ls <- lazy_dots(...)
  matrix_ls <- lazy_eval(ls[str_detect(names(ls), "x.")])
  names(matrix_ls) <- str_replace(names(matrix_ls), "x.", "")

    ni <- lapply(matrix_ls, function(matrix){
      nrow(matrix) - 1
    })

    p <- lapply(matrix_ls, function(matrix){
      ncol(matrix)
    })

    sample_covs <- lapply(matrix_ls, cov)

    scatterMat <- lapply(matrix_ls, A_func)

    overalln <- Reduce(`+`, ni)

    overallScatter <- Reduce(`+`, scatterMat) / overalln


    ahat2 <- (1 / ((overalln - 1) * (overalln + 2) * p[[1]])) * (tr(overallScatter %*% overallScatter) - (1 / overalln) * (tr(overallScatter)) ^ 2)

    npops <- length(matrix_ls)

    its <- combn(seq(npops), 2, simplify = FALSE)


    firstCSum <- Reduce(`+`, lapply(its, function(x){
      (ci_func(ni[[x[1]]], p[[1]]) + ci_func(ni[[x[2]]], p[[1]])) ^ 2
    }))

    secondCSum <- Reduce(`+`, lapply(ni, function(x){
      (ci_func(x, p[[1]])) ^ 2
    }))

    theta2 <- 4 * (ahat2 ^ 2) * (firstCSum + (npops - 1) * (npops - 2) * secondCSum)


  Schott2007_test.default(theta2, sample_covs)
}

#' @export
#'
Schott2007_test.default <- function(theta2, sample_covs){

  its <- combn(seq(length(sample_covs)), 2, simplify = FALSE)

  Reduce(`+`, lapply(its, function(x){
    (tr(sample_covs[[x[1]]] - sample_covs[[x[2]]]) ^ 4) / theta2
    }))

}
