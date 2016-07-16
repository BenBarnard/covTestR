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
dataDftoMatrix <- function(data){
  data %>%
    select(-Group) %>%
    acast(Subjects~Variables,
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
