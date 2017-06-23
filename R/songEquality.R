#' Song Conical Variates For Structure of Covariance Tests
#'
#' @param x
#' @param Sigma
#'
#' @return
#' @export
#'
#' @examples
songEquality <- function(x){
  UseMethod("songEquality")
}

#' @export
#' @keywords internal
songEquality.list <- function(x){
  sampleCov <- lapply(x, cov)
  covDiffs <- Reduce(cbind, lapply(sampleCov, function(x){x - sampleCov[[1]]})[-1])
  svdlist <- svd(covDiffs)
  scatters <- lapply(x, A_func)
  ns <- lapply(x, nrow)
  pooled <- overall_cov_func(scatters, ns)
  sumDiff <- Reduce(`+`, lapply(sampleCov, function(x){
    (x - pooled) %*% (x - pooled)
  }))

  weights <- apply(svdlist$u, 2, songHelper, Num = sumDiff, Dem = pooled)
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
