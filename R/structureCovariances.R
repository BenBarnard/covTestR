#' Structure of Covariances
#'
#' @param x data
#' @param Sigma Population covariance matrix
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
#' @examples structureCovariances(iris, group = Species)
structureCovariances <- function(x, Sigma, ..., covTest = Nagao1973){
  UseMethod("structureCovariances")
}

#' @export
#' @keywords internal
#' @importFrom dplyr as_data_frame
#' @importFrom stats setNames
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
structureCovariances.data.frame <- function(x, Sigma, ..., covTest = Nagao1973){
  dots <- lazy_dots(...)
  groupname <- names(unique(x[paste(dots$group$expr)]))
  group <- as.character(unique(x[[paste(dots$group$expr)]]))
  dots <- dots[!("group" %in% names(dots))]
  x <- setNames(lapply(group, function(y){
    as.matrix(x[x[groupname] == y,][names(x) != groupname])
  }), group)
  mat <- do.call(covTest, c(x = list(x), Sigma = list(Sigma), lazy_eval(dots)))
  mat
}

#' @export
#' @keywords internal
#' @importFrom dplyr as_data_frame
#' @importFrom stats setNames
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
structureCovariances.grouped_df <- function(x, Sigma, ..., covTest = Nagao1973){
  groups <- attributes(x)$labels
  x <- as_data_frame(x)
  group <- as.character(groups[,1])
  groupname <- names(groups)
  ls <- setNames(lapply(group, function(y){
    as.matrix(x[x[groupname] == y,][names(x) != groupname])
  }), group)
  dots <- lazy_dots(...)
  lapply(ls, function(x){
    mat <- do.call(covTest, c(x = list(x), Sigma = list(Sigma), lazy_eval(dots)))
    mat
  })
}

#' @export
#' @keywords internal
#' @importFrom dplyr as_data_frame
#' @importFrom stats setNames
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
structureCovariances.resample <- function(x, Sigma, ..., covTest = Nagao1973){
  x <- as_data_frame(x)
  groups <- attributes(x)$labels
  group <- as.character(groups[,1])
  groupname <- names(groups)
  ls <- setNames(lapply(group, function(y){
    as.matrix(x[x[groupname] == y,][names(x) != groupname])
  }), group)
  dots <- lazy_dots(...)
  lapply(ls, function(x){
    mat <- do.call(covTest, c(x = list(x), Sigma = list(Sigma), lazy_eval(dots)))
    mat
  })
}
