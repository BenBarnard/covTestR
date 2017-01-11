#' Turn no a Tidy data frame into data matrix (helper function)
#'
#' @param data not a tidy dataframe
#' @param group group
#' @param method method
#' @param ...  other options
#' @param .dots other options
#'
#' @importFrom plyr dlply
#' @importFrom lazyeval lazy_eval
#' @importFrom lazyeval lazy_dots
#'
#' @keywords internal
#'
#'
dataDftoMatrix <- function(data, group, method, ..., .dots){
  do.call(what = paste(method),
          args = c(x = dlply(.data = data,
                             .variables = group,
                             .fun = Tidy_,
                             group = group),
                   group = group,
                   lazy_eval(lazy_dots(...)),
                   lazy_eval(.dots)
          )
  )
}

#' not tidy helper
#'
#' @param data data
#' @param group grouping variable
#'
#' @importFrom dplyr select
#'
#' @keywords internal
#'
#'
Tidy_ <- function(data, group){
  as.matrix(select(data, -eval(group)))
}

#' Trace of Matrix
#'
#' @param mat matrix
#'
#' @keywords internal
#'
#'
tr <- function(mat){
  sum(diag(mat))
}
