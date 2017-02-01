#' Paste Wrapper
#'
#' @export
#' @keywords internal
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#'
past <- function(xmin, ..., xmax){
  x <- c(list(...))

  if(length(x[[1]]) > 0){
    namevec <- c(xmin, paste0(", ", x[[1]]), paste0(" and ", xmax))
  }else{
    namevec <- c(xmin, paste0(" and ", xmax))
  }
  namevec
}
