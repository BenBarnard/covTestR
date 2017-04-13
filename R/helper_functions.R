#' Use Method Creater for Tests
#' @keywords internal
#' @export
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#'
helper <- function(method){
  func <- function(x, ...){
    groups <- attributes(x)$labels
    x <- as.data.frame(x)
    dots <- lazy_dots(...)
    if(is.null(groups)){
      groupname <- names(unique(x[paste(dots$group$expr)]))
      group <- as.character(unique(x[paste(dots$group$expr)][,1]))
    }else{
      group <- as.character(groups[,1])
      groupname <- names(groups)
    }

    ls <- setNames(lapply(group, function(y){
      as.matrix(x[x[groupname] == y,][names(x) != groupname])
      }), group)

    do.call(what = method,
            args = list(x = ls,
                        .dots = dots))
  }
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
