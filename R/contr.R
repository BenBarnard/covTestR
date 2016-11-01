#' Title
#'
#' @param stat
#' @param control.var
#'
#' @return
#' @export
#'
#' @examples
contr <- function(stat, control.var){
  beta <- - cov(stat, control.var) / var(control.var)
  stat + beta * (control.var - 0)
}
