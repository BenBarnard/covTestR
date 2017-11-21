#' Test of Homogeneity of Covariance Matrices given by Ishii and Aoshima 2016
#'
#' @inherit homogeneityCovariances
#'
#' @references Ishii, A., Yata, K., and Aoshima, M. (2016). Asymptotic properties of the first pricipal
#' component and equality tests of covariance matrices in high-dimesion, low-sample-size context. Journal
#' of Statistical Planning and Inference, 170:186-199. \doi{10.1016/j.jspi.2015.10.007}
#'
#' @export
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#' @importFrom stats cov
#' @importFrom stats pf
#'
#' @examples 
#' irisSpecies <- unique(iris$Species)
#' 
#' iris_ls <- lapply(irisSpecies, 
#'     function(x){as.matrix(iris[iris$Species == x, 1:4])}
#'                  )
#'                  
#' names(iris_ls) <- irisSpecies
#' 
#' Ishii2016(iris_ls)
Ishii2016 <- function(x, ...) {

  ls <- lazy_dots(...)
  matrix_ls <- x

  statistic <- Ishii2016Stat(matrix_ls)

  xmin <- names(matrix_ls[1])
  xmax <- names(matrix_ls[length(matrix_ls)])
  xother <- names(matrix_ls[-c(1, length(matrix_ls))])

  data.name <- Reduce(paste0, past(xmin = xmin, xother, xmax = xmax))

  names(statistic) <- "F"

  parameter <- c(length(matrix_ls), 1)
  names(parameter) <- c("df1", "df2")

  null.value <- 0
  names(null.value) <- "difference in covariance matrices"

  p.value <- 1 - pf(statistic, df1 = parameter[1], df2 = parameter[2])


  obj <- list(statistic = statistic,
              parameter = parameter,
              p.value = p.value,
              estimate = NULL,
              null.value = null.value,
              alternative = "two.sided",
              method = "Ishii and Aoshima 2016 Homogeneity of Covariance Matrices Test",
              data.name = data.name)
  class(obj) <- "htest"
  obj
  }
