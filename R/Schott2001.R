#' Test of Equality of Covariances given by Schott 2001
#'
#' Performs 2 and k sample equality of covariance matrix test using Schott 2001
#'
#' @inheritParams Chaipitak2013
#'
#' @return Test statistic of the hypothesis test
#'
#' @export
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#' @importFrom stats cov
#' @importFrom stats pchisq
#'
#' @examples Schott2001(iris, group = Species)
#'
Schott2001 <- function(x, ...) {

  ls <- lazy_dots(...)
  matrix_ls <- x

  if(!("covariance" %in% class(x[[1]])) & ("matrix" %in% class(x[[1]]))){
    ns <- lapply(matrix_ls, function(matrix){
      nrow(matrix)
    })

    A_ls <- lapply(matrix_ls, A_func)

    sample_covs <- lapply(matrix_ls, cov)
  }

  if("covariance" %in% class(x[[1]])){
    ns <- lapply(matrix_ls, function(matrix){
      attributes(matrix)$df + 1
    })

    sample_covs <- matrix_ls

    A_ls <- mapply(function(sample_covs, ns){
      sample_covs * (ns - 1)
    }, sample_covs = sample_covs, ns = ns, SIMPLIFY = FALSE)
  }

  groups <- 1:length(matrix_ls)

  overall_cov <- overall_cov_func(A_ls, ns)
  invOverall_cov <- solve(overall_cov)

  overall_n <- Reduce(`+`, ns)

  xmin <- names(matrix_ls[1])
  xmax <- names(matrix_ls[length(matrix_ls)])
  xother <- names(matrix_ls[-c(1, length(matrix_ls))])

  data.name <- Reduce(paste0, past(xmin = xmin, xother, xmax = xmax))

  statistic <- (overall_n / 2) * (Reduce(`+`, lapply(groups, function(x){
    (ns[[x]] / overall_n) *
      tr(sample_covs[[x]] %*% invOverall_cov %*%
           sample_covs[[x]] %*% invOverall_cov)
  })) -
    Reduce(`+`, lapply(groups, function(z){
      Reduce(`+`, lapply(groups, function(y){
        (ns[[z]] * ns[[y]] / (overall_n ^ 2)) *
          tr(sample_covs[[z]] %*% invOverall_cov %*%
               sample_covs[[y]] %*% invOverall_cov)
      }))
    })))




  names(statistic) <- "Chi Squared"

  parameter <- 1
  names(parameter) <- "df"

  null.value <- 0
  names(null.value) <- "difference in covariances"

  p.value <- 1 - pchisq(statistic, parameter)

  estimate <- sample_covs
  names(estimate) <- paste0("covariance of ", names(matrix_ls))

  estimate <- if(nrow(estimate[[1]]) > 5){
    NULL
  }else{
    estimate
  }

  obj <- list(statistic = statistic,
              parameter = parameter,
              p.value = p.value,
              estimate = estimate,
              null.value = null.value,
              alternative = "two.sided",
              method = "Schott 2007 Equality of Covariance Test",
              data.name = data.name)
  class(obj) <- "htest"
  obj
}
