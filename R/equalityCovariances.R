#' Equality of Covariances
#'
#' @param x data
#' @param ... other options passed to covTest method
#' @param covTest equality of covariances test method
#'
#' @return A list with class "htest" containing the following components:
#'
#'\tabular{ll}{
#' \code{statistic} \tab the value of equality of covariance test statistic \cr
#' \tab \cr
#' \code{parameter} \tab the degrees of freedom for the chi-squared statistic \cr
#' \tab \cr
#' \code{p.value} \tab the p=value for the test \cr
#' \tab \cr
#' \code{estimate} \tab the estimated covariances if less than 5 dimensions \cr
#' \tab \cr
#' \code{null.value} \tab the specified hypothesized value of the covariance difference \cr
#' \tab \cr
#' \code{alternative} \tab a character string describing the alternative hyposthesis \cr
#' \tab \cr
#' \code{method} \tab a character string indicating what type of equality of covariance test was performed \cr
#' \tab \cr
#' \code{data.name} \tab a character string giving the names of the data
#'}
#'
#' @export
#'
#' @examples equalityCovariances(iris, group = Species)
equalityCovariances <- function(x, ..., covTest = BoxesM){
  UseMethod("equalityCovariances")
}

#' @export
#' @keywords internal
#' @importFrom dplyr as_data_frame
#' @importFrom stats setNames
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
equalityCovariances.data.frame <- function(x, ..., covTest = BoxesM){
  dots <- lazy_dots(...)
  groupname <- names(unique(x[paste(dots$group$expr)]))
  group <- as.character(unique(x[[paste(dots$group$expr)]]))
  dots <- dots[!("group" %in% names(dots))]
  x <- setNames(lapply(group, function(y){
    as.matrix(x[x[groupname] == y,][names(x) != groupname])
  }), group)
   mat <- do.call(covTest, c(x = list(x), lazy_eval(dots)))
   mat
}

#' @export
#' @keywords internal
#' @importFrom dplyr as_data_frame
#' @importFrom stats setNames
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
equalityCovariances.grouped_df <- function(x, ..., covTest = BoxesM){
  groups <- attributes(x)$labels
  x <- as_data_frame(x)
  group <- as.character(groups[,1])
  groupname <- names(groups)
  x <- setNames(lapply(group, function(y){
    as.matrix(x[x[groupname] == y,][names(x) != groupname])
  }), group)
  dots <- lazy_dots(...)
  mat <- do.call(covTest, c(x = list(x), lazy_eval(dots)))
  mat
}

#' @export
#' @keywords internal
#' @importFrom dplyr as_data_frame
#' @importFrom stats setNames
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
equalityCovariances.resample <- function(x, ..., covTest = BoxesM){
  x <- as_data_frame(x)
  groups <- attributes(x)$labels
  group <- as.character(groups[,1])
  groupname <- names(groups)
  ls <- setNames(lapply(group, function(y){
    as.matrix(x[x[groupname] == y,][names(x) != groupname])
  }), group)
  dots <- lazy_dots(...)
  lapply(ls, function(x){
    mat <- do.call(covTest, c(x = list(x), lazy_eval(dots)))
    mat
  })
}

#' @export
#' @keywords internal
#' @importFrom dplyr as_data_frame
#' @importFrom stats setNames
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
equalityCovariances.list <- function(x, ..., covTest = BoxesM){
  ls <- lazy_dots(...)
  matrix_ls <- x
    mat <- do.call(covTest, c(x = matrix_ls, lazy_eval(dots)))
    mat
}
