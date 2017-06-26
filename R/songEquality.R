#' Song Conical Variates For Structure of Covariance Tests
#'
#' @param x data
#' @param Sigma covariance matrix
#'
#' @export
#' @keywords internal
songEquality <- function(x){
  UseMethod("songEquality")
}

#' @export
#' @keywords internal
songEquality.list <- function(x){
  sampleCov <- lapply(x, cov)
  covDiffs <- Reduce(cbind, lapply(sampleCov, function(x){x - sampleCov[[1]]})[-1])
  svdlist <- svd(covDiffs)
  ns <- lapply(x, nrow)
  scatters <- mapply(`*`, sampleCov, lapply(ns, function(x){x - 1}), SIMPLIFY = FALSE)
  pooled <- Reduce(`+`, scatters) / Reduce(`+`, lapply(ns, function(x){x - 1}))
  sumDiff <- Reduce(`+`, lapply(sampleCov, function(x){
    (x - pooled) %*% (x - pooled)
  }))
  weights <- apply(svdlist$u, 2, songHelper, Num = sumDiff, Dem = pooled)
  ordering <- order(weights, decreasing = TRUE)
  d <- weights[ordering]
  u <- svdlist$u[, ordering]
  list(u = u, d = d)
}

#' Song Helper
#' @export
#' @keywords internal
songHelper <- function(ai, Num, Dem){
  (t(ai) %*% Num %*% ai) / (t(ai) %*% Dem %*% ai)
}
