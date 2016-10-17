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
#' @examples
#'
Ishii2016_test <- function(data, ...) {
  UseMethod("Ishii2016_test")
}

#' @export
#'
#' @importFrom lazyeval expr_find
#'
Ishii2016_test.data.frame <- function(x, group, ..., variables, samples, value, tidy = FALSE){
  if(tidy == TRUE){
    tidyDataDftoMatrix(data = x,
                       group = expr_find(group),
                       variables = expr_find(variable),
                       samples = expr_find(samples),
                       value = expr_find(value),
                       test = expr_find(Ishii2016_test.matrix))
  }else{
    dataDftoMatrix(data = x,
                   group = expr_find(group),
                   test = expr_find(Ishii2016_test.matrix))
  }
}

#' @export
#'
#' @importFrom plyr llply
#' @importFrom plyr mlply
#'
Ishii2016_test.matrix<- function(...){

}

#' @export
#'
Ishii2016_test.default <- function(lambdaTilde, gammaTilde, hTilde){
lambdaTilde1 <- lambdaTilde[[1]]
lambdaTilde2 <- lambdaTilde[[2]]

(lambdaTilde1 / lambdaTilde2) * hTilde * gammaTilde
}
