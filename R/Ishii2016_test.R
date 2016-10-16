#' Test of Equality of Covariances given by Schott 2007
#'
#' @param data tidy data frame
#' @param ... other
#'
#'
#'
#' @return Test Statistic for Schott 2007
#' @export
#'
#' @examples
#'
Ishii2016_test <- function(data, ...) {
  UseMethod("Ishii2016_test")
}



#' @export
#'
#' @importFrom plyr dlply
#' @importFrom plyr llply
#' @importFrom plyr mlply
#' @importFrom plyr .
#' @importFrom magrittr %>%
#'
Ishii2016_test.matrix<- function(...){

}



#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom plyr dlply
#' @importFrom plyr .
#'
Ishii2016_test.data.frame <- function(data, ...){

}


#' @export
#'
Ishii2016_test.default <- function(n, p, ahat2, ahat2i, sample_covs){

}
