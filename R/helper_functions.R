#' Turn Tidy data frame into data matrix (helper function)
#'
#' @param data tidy dataframe
#'
#' @importFrom plyr dlply
#'
#' @keywords internal
#'
#' @export
#'
tidyDataDftoMatrix <- function(data, group, others, test){
  browser()
  do.call(what = paste(test),
          args = dlply(.data = data,
                       .variables = group,
                       .fun = tidy_,
                       group = group,
                       others = others
          )
  )
}

#' Tidy helper
#'
#' @param data
#' @param group
#' @param others
#' @param test
#'
#' @importFrom reshape2 acast
#' @importFrom dplyr select
#'
#' @return
#' @export
#'
#' @examples
tidy_ <- function(data, group, others){
  as.matrix(acast(select(data,
               -eval(group)),
        eval(others$samples$expr)~eval(others$variable$expr),
        value.var = paste(others$value$expr)))
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
