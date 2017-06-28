sdiff <- function(ls){
  covs <- lapply(ls, cov)
  diffs <- lapply(covs, function(x){x - covs[[1]]})[-1]
  redMat <- svd(Reduce(cbind, diffs))$u
  redMat
}

sdiffpool <- function(ls){
  covs <- lapply(ls, cov)
  scatters <- mapply(`*`, covs, lapply(ls, function(x){nrow(x) - 1}), SIMPLIFY = FALSE)
  pooled <- Reduce(`+`, scatters) / Reduce(`+`, lapply(ls, function(x){nrow(x) - 1}))
  diffs <- lapply(covs, function(x){x - pooled})
  redMat <- svd(Reduce(cbind, diffs))$u
  redMat
}

sconcat <- function(ls){
  covs <- lapply(ls, cov)
  redMat <- svd(Reduce(cbind, covs))$u
  redMat
}

reorderSdiff <- function(ls){
  sampleCov <- lapply(ls, cov)
  covDiffs <- Reduce(cbind, lapply(sampleCov, function(x){x - sampleCov[[1]]})[-1])
  svdlist <- svd(covDiffs)
  ns <- lapply(ls, nrow)
  scatters <- mapply(`*`, sampleCov, lapply(ns, function(x){x - 1}), SIMPLIFY = FALSE)
  pooled <- Reduce(`+`, scatters) / Reduce(`+`, lapply(ns, function(x){x - 1}))
  sumDiff <- Reduce(`+`, lapply(sampleCov, function(x){
    (x - pooled) %*% (x - pooled)
  }))

  weights <- apply(svdlist$u, 2, reorderHelper, Num = sumDiff, Dem = pooled)
  ordering <- order(weights, decreasing = TRUE)
  d <- weights[ordering]
  u <- svdlist$u[, ordering]
  u
}




reorderHelper <- function(ai, Num, Dem){
  (t(ai) %*% Num %*% ai) / (t(ai) %*% Dem %*% ai)
}
