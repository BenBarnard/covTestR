#' Turn no a Tidy data frame into data matrix (helper function)
#'
#' @param data not a tidy dataframe
#'
#' @importFrom plyr dlply
#'
#' @keywords internal
#'
#' @export
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
#' @export
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
#' @export
#'
tr <- function(mat){
  sum(diag(mat))
}
