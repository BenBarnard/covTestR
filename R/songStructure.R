#' Song Conical Variates For Structure of Covariance Tests
#'
#' @param x
#' @param Sigma
#'
#' @return
#' @export
#'
#' @examples
songStructure <- function(x, Sigma){
  UseMethod("songStructure")
}


#' @export
#' @keywords internal
songStructure.covariance <- function(x, Sigma){
  svdlist <- svd(x)
  weights <- apply(svdlist$u, 2, songHelper, Num = (sampleCov - Sigma) %*% (sampleCov - Sigma), Dem = Sigma)
  ordering <- order(weights, decreasing = TRUE)
  d <- weights[ordering]
  u <- svdlist$u[, ordering]
  list(u = u, d = d)
}

#' @export
#' @keywords internal
songStructure.matrix <- function(x, Sigma){
  sampleCov <- cov(x)
  svdlist <- svd(sampleCov)
  weights <- apply(svdlist$u, 2, songHelper, Num = (sampleCov - Sigma) %*% (sampleCov - Sigma), Dem = Sigma)
  ordering <- order(weights, decreasing = TRUE)
  d <- weights[ordering]
  u <- svdlist$u[, ordering]
  list(u = u, d = d)
}

#' @export
#' @keywords internal
songHelper <- function(ai, Num, Dem){
  (t(ai) %*% Num %*% ai) / (t(ai) %*% Dem %*% ai)
}
