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
dataDftoMatrix <- function(data, group, test){
  do.call(what = paste(test),
          args = c(x = dlply(.data = data,
                       .variables = group,
                       .fun = Tidy_,
                       group = group),
                   group = group)
  )
}

#' tidy helper
#'
#' @param data
#' @param group
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
