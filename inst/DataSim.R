library(MASS)
library(pushoverr)
set_pushover_user(user = "ufmfa6vc9s2fc2eh6phop9ej5ebxum")
set_pushover_app(token = "azrd3hwrwgh2gs6igbvb8yy4mftoi7")


## Autoregressive

maxdimensions <- 500
maxSampleSize <- 100
replications <- 1000

sigma1 <- sigma3 <- matrix(0, nrow = maxdimensions, ncol = maxdimensions)
for(i in 1:maxdimensions){
  for(j in 1:maxdimensions){
    sigma1[i, j] <- 0.1 ^ (abs(i - j))
  }
}

for(i in 1:maxdimensions){
  for(j in 1:maxdimensions){
    sigma3[i, j] <- 0.3 ^ (abs(i - j))
  }
}

Sigmaj <- list("Zero" = sigma1,
               "One" = sigma3,
               "Two" = sigma1)

save(Sigmaj, file = "E:/Ben/Box Sync/Statistics/Dissertation/mvnData/Autoregressive/Sigmaj.RData")

mvndata <- replicate(replications,
                     list("Zero1" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                            Sigma = Sigmaj[names(Sigmaj) == "Zero"][[1]][1:maxdimensions, 1:maxdimensions]),
                          "Zero2" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                            Sigma = Sigmaj[names(Sigmaj) == "Zero"][[1]][1:maxdimensions, 1:maxdimensions]),
                          "Zero3" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                            Sigma = Sigmaj[names(Sigmaj) == "Zero"][[1]][1:maxdimensions, 1:maxdimensions]),
                          "One" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                          Sigma = Sigmaj[names(Sigmaj) == "One"][[1]][1:maxdimensions, 1:maxdimensions]),
                          "Two" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                          Sigma = Sigmaj[names(Sigmaj) == "Two"][[1]][1:maxdimensions, 1:maxdimensions])),
                     simplify = FALSE)

save(mvndata, file = "E:/Ben/Box Sync/Statistics/Dissertation/mvnData/Autoregressive/mvndata.RData")

pushover(message = "Autoregressive",
         title = "dataSim")

## Compound Symmetry

Sigmaj <- list("Zero" = .99 * diag(1, maxdimensions) + 0.01 * rep(1, maxdimensions) %*% t(rep(1, maxdimensions)),
               "One" = .95 * diag(1, maxdimensions) + 0.05 * rep(1, maxdimensions) %*% t(rep(1, maxdimensions)),
               "Two" = .99 * diag(1, maxdimensions) + 0.01 * rep(1, maxdimensions) %*% t(rep(1, maxdimensions)))

save(Sigmaj, file = "E:/Ben/Box Sync/Statistics/Dissertation/mvnData/CompoundSymmetry/Sigmaj.RData")

mvndata <- replicate(replications,
                     list("Zero1" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                            Sigma = Sigmaj[names(Sigmaj) == "Zero"][[1]][1:maxdimensions, 1:maxdimensions]),
                          "Zero2" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                            Sigma = Sigmaj[names(Sigmaj) == "Zero"][[1]][1:maxdimensions, 1:maxdimensions]),
                          "Zero3" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                            Sigma = Sigmaj[names(Sigmaj) == "Zero"][[1]][1:maxdimensions, 1:maxdimensions]),
                          "One" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                          Sigma = Sigmaj[names(Sigmaj) == "One"][[1]][1:maxdimensions, 1:maxdimensions]),
                          "Two" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                          Sigma = Sigmaj[names(Sigmaj) == "Two"][[1]][1:maxdimensions, 1:maxdimensions])),
                     simplify = FALSE)

save(mvndata, file = "E:/Ben/Box Sync/Statistics/Dissertation/mvnData/CompoundSymmetry/mvndata.RData")

pushover(message = "Compound Symmetry",
         title = "dataSim")


## Identity

Sigmaj <- list("Zero" = diag(rep(1, maxdimensions)),
               "One" = diag(rep(1.5, maxdimensions)),
               "Two" = diag(rep(1, maxdimensions)))

mvndata <- replicate(replications,
                     list("Zero1" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                            Sigma = Sigmaj[names(Sigmaj) == "Zero"][[1]][1:maxdimensions, 1:maxdimensions]),
                          "Zero2" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                            Sigma = Sigmaj[names(Sigmaj) == "Zero"][[1]][1:maxdimensions, 1:maxdimensions]),
                          "Zero3" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                            Sigma = Sigmaj[names(Sigmaj) == "Zero"][[1]][1:maxdimensions, 1:maxdimensions]),
                          "One" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                          Sigma = Sigmaj[names(Sigmaj) == "One"][[1]][1:maxdimensions, 1:maxdimensions]),
                          "Two" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                          Sigma = Sigmaj[names(Sigmaj) == "Two"][[1]][1:maxdimensions, 1:maxdimensions])),
                     simplify = FALSE)

save(mvndata, file = "E:/Ben/Box Sync/Statistics/Dissertation/mvnData/Identity/mvndata.RData")

pushover(message = "Identity",
         title = "dataSim")



## Sri2014

sigma <- diag(1 + (-1) ^ (seq(1, maxdimensions) + 1) * runif(maxdimensions, 0, 1))
rhomat1 <- rhomat3 <- matrix(0, nrow = maxdimensions, ncol = maxdimensions)

for(i in 1:maxdimensions){
  for(j in 1:maxdimensions){
    rhomat1[i, j] <- 0.1 ^ (abs(i - j) ^ (1 / 10))
  }
}

for(i in 1:maxdimensions){
  for(j in 1:maxdimensions){
    rhomat3[i, j] <- 0.3 ^ (abs(i - j) ^ (1 / 10))
  }
}

Sigmaj <- list("Zero" = sigma %*% rhomat1 %*% sigma,
               "One" = sigma %*% rhomat3 %*% sigma,
               "Two" = sigma %*% rhomat1 %*% sigma)

save(Sigmaj, file = "E:/Ben/Box Sync/Statistics/Dissertation/mvnData/Sri2014/Sigmaj.RData")

mvndata <- replicate(replications,
                     list("Zero1" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                            Sigma = Sigmaj[names(Sigmaj) == "Zero"][[1]][1:maxdimensions, 1:maxdimensions]),
                          "Zero2" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                            Sigma = Sigmaj[names(Sigmaj) == "Zero"][[1]][1:maxdimensions, 1:maxdimensions]),
                          "Zero3" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                            Sigma = Sigmaj[names(Sigmaj) == "Zero"][[1]][1:maxdimensions, 1:maxdimensions]),
                          "One" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                          Sigma = Sigmaj[names(Sigmaj) == "One"][[1]][1:maxdimensions, 1:maxdimensions]),
                          "Two" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                          Sigma = Sigmaj[names(Sigmaj) == "Two"][[1]][1:maxdimensions, 1:maxdimensions])),
                     simplify = FALSE)

save(mvndata, file = "E:/Ben/Box Sync/Statistics/Dissertation/mvnData/Sri2014/mvndata.RData")

pushover(message = "Sri2014",
         title = "dataSim")



## Unstructured

omega <- runif(max(dimensions), 1, 5)
j <- setNames(c(0, 1, 2), c("Zero", "One", "Two"))

deltaj <- lapply(j, function(j){
  deltaj <- matrix(0, max(dimensions), max(dimensions))
  for(a in 1:max(dimensions)){
    for(b in 1:max(dimensions)){
      deltaj[a, b] <- ((-1) ^ (a + b)) * ((0.2 * (j + 2)) ^ (abs(a - b) ^ (1 / 10)))
    }
  }
  deltaj
})

Sigmaj <- lapply(deltaj, function(deltaj){
  diag(omega) %*% deltaj %*% diag(omega)
})

save(Sigmaj, file = "E:/Ben/Box Sync/Statistics/Dissertation/mvnData/Unstructured/Sigmaj.RData")

mvndata <- replicate(replications,
                     list("Zero1" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                            Sigma = Sigmaj[names(Sigmaj) == "Zero"][[1]][1:maxdimensions, 1:maxdimensions]),
                          "Zero2" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                            Sigma = Sigmaj[names(Sigmaj) == "Zero"][[1]][1:maxdimensions, 1:maxdimensions]),
                          "Zero3" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                            Sigma = Sigmaj[names(Sigmaj) == "Zero"][[1]][1:maxdimensions, 1:maxdimensions]),
                          "One" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                          Sigma = Sigmaj[names(Sigmaj) == "One"][[1]][1:maxdimensions, 1:maxdimensions]),
                          "Two" = mvrnorm(n = maxSampleSize, mu = rep(0, maxdimensions),
                                          Sigma = Sigmaj[names(Sigmaj) == "Two"][[1]][1:maxdimensions, 1:maxdimensions])),
                     simplify = FALSE)

save(mvndata, file = "E:/Ben/Box Sync/Statistics/Dissertation/mvnData/Unstructured/mvndata.RData")

pushover(message = "Unstructured",
         title = "dataSim")
