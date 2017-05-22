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
songStructure.covariance <- function(x){
  svdlist <- svd(x)
  di <- diag(1, nrow(sampleCov))
  weights <- apply(svdlist$u, 2, songHelper, Num = (sampleCov - di) %*% (sampleCov - di))
  ordering <- order(weights, decreasing = TRUE)
  d <- weights[ordering]
  u <- svdlist$u[, ordering]
  list(u = u, d = d)
}

#' @export
#' @keywords internal
songStructure.matrix <- function(x){
  sampleCov <- cov(x)
  svdlist <- svd(sampleCov)
  di <- diag(1, nrow(sampleCov))
  weights <- apply(svdlist$u, 2, songStructureHelper, Num = (sampleCov - di) %*% (sampleCov - di))
  ordering <- order(weights, decreasing = TRUE)
  d <- weights[ordering]
  u <- svdlist$u[, ordering]
  list(u = u, d = d)
}

#' @export
#' @keywords internal
songStructureHelper <- function(ai, Num){
  t(ai) %*% Num %*% ai
}
