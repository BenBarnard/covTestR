#' Frobenius Norm
#'
#' @param x data
#'
#' @export
#'
#' @keywords internal
frobenius_norm <- function(x){
  x$`Frobenius Norm` <- 0
x[x$type == "identity" & x$populations == 3,]$`Frobenius Norm` <-
  ddply(x[x$type == "identity" & x$populations == 3,], .(ReductionMethod, originaldimensions, ReducedDimension, test, difference, Samples),
        function(x){
          2 * norm(diag(1,x$originaldimensions) -
                     diag(1 +x$difference, x$originaldimensions)) / x$originaldimensions
        })$V1
x[x$type == "toeplitz" & x$populations == 3,]$`Frobenius Norm` <-
  ddply(x[x$type == "toeplitz" & x$populations == 3,], .(ReductionMethod, originaldimensions, ReducedDimension, test, difference, Samples),
        function(x){
          2 *
            norm(toeplitz(.5 ^ seq(0, (x$originaldimensions - 1))) -
                   (toeplitz(.5 ^ seq(0, (x$originaldimensions - 1))) + x$difference)) / x$originaldimensions
        })$V1
x[x$type == "identity" & x$populations == 2,]$`Frobenius Norm` <-
  ddply(x[x$type == "identity" & x$populations == 2,], .(ReductionMethod, originaldimensions, ReducedDimension, test, difference, Samples),
        function(x){
          norm(diag(1,x$originaldimensions) -
                 diag(1 +x$difference, x$originaldimensions)) / x$originaldimensions
        })$V1
x[x$type == "toeplitz" & x$populations == 2,]$`Frobenius Norm` <-
  ddply(x[x$type == "toeplitz" & x$populations == 2,], .(ReductionMethod, originaldimensions, ReducedDimension, test, difference, Samples),
        function(x){
          norm(toeplitz(.5 ^ seq(0, (x$originaldimensions - 1))) -
                 (toeplitz(.5 ^ seq(0, (x$originaldimensions - 1))) + x$difference)) / x$originaldimensions
        })$V1
x
   }



#' Title
#'
#' @param x data
#'
#' @export
#'
#' @keywords internal
frobenius_norm2 <- function(x){
  x$`Frobenius Norm` <- 0
  if(any(x$type == "elliptical")){
    x[x$type == "elliptical" & x$populations == 3,]$`Frobenius Norm` <-
    ddply(x[x$type == "elliptical" & x$populations == 3,], .(ReductionMethod, originaldimensions, ReducedDimension, test, difference, Samples),
          function(x){
            2 * norm(cov_maker(keepers = list(c(20, 1, rep(5, 8))),
                                                            offs = list(0, nrow = 10, ncol = x$originaldimensions - 10),
                                                             losers = list(c(1, rep(0, x$originaldimensions - 11)))) -
                                                     (cov_maker(keepers = list(c(20, 1, rep(5, 8))),
                                                                offs = list(0, nrow = 10, ncol = x$originaldimensions - 10),
                                                                losers = list(c(1, rep(0, x$originaldimensions - 11)))) + x$difference)) / x$originaldimensions
          })$V1
  x[x$type == "elliptical" & x$populations == 2,]$`Frobenius Norm` <-
    ddply(x[x$type == "elliptical" & x$populations == 2,], .(ReductionMethod, originaldimensions, ReducedDimension, test, difference, Samples),
          function(x){
            norm(cov_maker(keepers = list(c(20, 1, rep(5, 8))),
                           offs = list(0, nrow = 10, ncol = x$originaldimensions - 10),
                           losers = list(c(1, rep(0, x$originaldimensions - 11)))) -
                   (cov_maker(keepers = list(c(20, 1, rep(5, 8))),
                              offs = list(0, nrow = 10, ncol = x$originaldimensions - 10),
                              losers = list(c(1, rep(0, x$originaldimensions - 11)))) + x$difference)) / x$originaldimensions
          })$V1
  }

  if(any(x$type == "toeplitz")){
  x[x$type == "toeplitz" & x$populations == 3,]$`Frobenius Norm` <-
    ddply(x[x$type == "toeplitz" & x$populations == 3,], .(ReductionMethod, originaldimensions, ReducedDimension, test, difference, Samples),
          function(x){
            2 *
              norm(toeplitz(.5 ^ seq(0, (x$originaldimensions - 1))) -
                     (toeplitz(.5 ^ seq(0, (x$originaldimensions - 1))) * x$difference)) / x$originaldimensions
          })$V1

  x[x$type == "toeplitz" & x$populations == 2,]$`Frobenius Norm` <-
    ddply(x[x$type == "toeplitz" & x$populations == 2,], .(ReductionMethod, originaldimensions, ReducedDimension, test, difference, Samples),
          function(x){
            norm(toeplitz(.5 ^ seq(0, (x$originaldimensions - 1))) -
                   (toeplitz(.5 ^ seq(0, (x$originaldimensions - 1))) * x$difference)) / x$originaldimensions
          })$V1
  }
  x
}
