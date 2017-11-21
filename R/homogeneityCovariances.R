#' Homogeneity of Covariance Matrices Test
#' 
#' Performs 2 and k sample homogeneity of covariance matrices test using test, 
#' 'covTest.'
#' 
#' @param x data as a data frame, list of matrices, grouped data frame, or resample object
#' @param ... other options passed to covTest method
#' @param covTest homogeneity of covariance matrices test method
#'
#' @return A list with class "htest" containing the following components:
#'
#'\tabular{ll}{
#' \code{statistic} \tab the value of homogeneity of covariance test statistic \cr
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
#' \code{method} \tab a character string indicating what type of homogeneity of covariance test was performed \cr
#' \tab \cr
#' \code{data.name} \tab a character string giving the names of the data
#'}
#'
#' @details The \code{\link{homogeneityCovariances}} function is a wrapper function that formats the data 
#'   for the specific \code{covTest} functions.
#'
#' @export
#'
#' @examples homogeneityCovariances(iris, group = Species)
homogeneityCovariances <- function(x, ..., covTest = BoxesM){
  UseMethod("homogeneityCovariances")
}

#' @export
#' @keywords internal
#' @importFrom dplyr as_data_frame
#' @importFrom stats setNames
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
homogeneityCovariances.data.frame <- function(x, ..., covTest = BoxesM){
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
homogeneityCovariances.grouped_df <- function(x, ..., covTest = BoxesM){
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
homogeneityCovariances.resample <- function(x, ..., covTest = BoxesM){
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
homogeneityCovariances.list <- function(x, ..., covTest = BoxesM){
  dots <- lazy_dots(...)
  matrix_ls <- x
    mat <- do.call(covTest, c(x = matrix_ls, lazy_eval(dots)))
    mat
}
