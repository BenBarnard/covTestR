#' Title
#'
#' @param stat
#' @param control.var
#'
#' @return
#' @export
#'
#' @examples
contr <- function(stat, control.var, quantile_variate){
  df <- data.frame(stat = stat, control.var = control.var)
  quant <- quantile(df$control.var, probs <- quantile_variate)

  browser()

  n <- length(df$stat)

  empcdf <- sapply(df$stat, function(y){
    nrow(df[df$stat <= y,]) / n
  })

  N01 <- sapply(df$stat, function(y){
    nrow(df[df$control.var <= quant[[1]] & df$stat > y,])
  })

  N00 <- sapply(df$stat, function(y){
    nrow(df[df$control.var <= quant[[1]] & df$stat <= y,])
    })

  N10 <- sapply(df$stat, function(y){
    nrow(df[df$control.var > quant[[1]] & df$stat <= y,])
  })

  N11 <- sapply(df$stat, function(y){
    nrow(df[df$control.var > quant[[1]] & df$stat > y,])
  })

  N0 <- nrow(df[df$control.var <= quant[[1]],])
  N1 <- nrow(df[df$control.var > quant[[1]],])

  p00 <- (quantile_variate * N00) / (N0)
  p10 <- ((1 - quantile_variate) * N10) / (N1)
  cdf <- p00 + p10
  cdf
}
