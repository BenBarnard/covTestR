#' Turn Tidy data frame into data matrix (helper function)
#'
#' @param data tidy dataframe
#'
#' @importFrom magrittr %>%
#' @importFrom reshape2 acast
#' @importFrom dplyr select
#'
#' @keywords internal
#'
#' @export
#'
tidyDataDftoMatrix <- function(data, group, variables, samples){
  acast(select(data,
               -expr_find(group)),
        expr_find(samples)~expr_find(variables),
        value.var = "Value")
}

#' Turn Tidy data frame into data matrix (helper function)
#'
#' @param data tidy dataframe
#'
#' @importFrom magrittr %>%
#' @importFrom reshape2 acast
#' @importFrom dplyr select
#'
#' @keywords internal
#'
#' @export
#'
dataDftoMatrix <- function(data, group){
  acast(select(data,
               -expr_find(group)),
        Subjects~Variables,
        value.var = "Value")
}

#' Trace of Matrix
#'
#' @param mat matrix
#'
#' @keywords internal
#'
#' @export
#'
tr <- function(mat){
  sum(diag(mat))
}


function(x, test, group, variables, samples){
  if(tidy == TRUE){
  do.call(what = Chaipitak2013_test.matrix,
          args = c(dlply(.data = x,
                         .variables = expr_find(group),
                         .fun = dataDftoMatrix),
                   group = expr_find(group),
                   variables = expr_find(variables),
                   samples = expr_find(samples)))
}else{
  browser()

}
}
