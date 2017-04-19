#' Use Method Creater for Tests
#' @keywords internal
#' @export
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stats setNames
#' @importFrom dplyr as_data_frame
#'
helper <- function(method){
  func <- function(x, ...){
    if(is.null(attributes(x)$labels)){
      x <- dplyr::as_data_frame(x)
      groups <- attributes(x)$labels
    }else{
      groups <- attributes(x)$labels
      x <- dplyr::as_data_frame(x)
    }
    dots <- lazyeval::lazy_dots(...)
    if(is.null(groups)){
      groupname <- names(unique(x[paste(dots$group$expr)]))
      group <- as.character(unique(x[[paste(dots$group$expr)]]))
    }else{
      group <- as.character(groups[,1])
      groupname <- names(groups)
    }

    ls <- stats::setNames(lapply(group, function(y){
      as.matrix(x[x[groupname] == y,][names(x) != groupname])
    }), group)

    do.call(what = method,
            args = list(x = ls,
                        .dots = dots))
  }
}

#' Use Method Creater for Tests
#' @keywords internal
#' @export
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stats setNames
#' @importFrom dplyr as_data_frame
#'
helperOne <- function(method){
  func <- function(x, ...){
    ls <- as.matrix(x)
    dots <- lazyeval::lazy_dots(...)
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
