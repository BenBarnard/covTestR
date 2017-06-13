library(MASS)
library(plyr)
library(EqualCov)
library(dplyr)
library(pushoverr)
set_pushover_user(user = "ufmfa6vc9s2fc2eh6phop9ej5ebxum")
set_pushover_app(token = "azrd3hwrwgh2gs6igbvb8yy4mftoi7")

dimensions <- c(20, 40, 80, 160)
SampleSize <- c(5, 10, 20, 40)
replications <- 1000
grid <- expand.grid(dimensions = dimensions, SampleSize = SampleSize, Sigmaj = c("Zero", "One", "Two"))


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

save(Sigmaj, file = "E:/Ben/Box Sync/Statistics/Unstructured/Sigmaj.RData")

mvndata <- mapply(function(SampleSize, Sig, Sigmaj, dimensions, replications){
  df <- replicate(replications,
                  list("Zero1" = mvrnorm(n = SampleSize, mu = rep(0, dimensions),
                                         Sigma = Sigmaj[names(Sigmaj) == "Zero"][[1]][1:dimensions, 1:dimensions]),
                       "Zero2" = mvrnorm(n = SampleSize, mu = rep(0, dimensions),
                                         Sigma = Sigmaj[names(Sigmaj) == "Zero"][[1]][1:dimensions, 1:dimensions]),
                       "Zero3" = mvrnorm(n = SampleSize, mu = rep(0, dimensions),
                                         Sigma = Sigmaj[names(Sigmaj) == "Zero"][[1]][1:dimensions, 1:dimensions]),
                       "One" = mvrnorm(n = SampleSize, mu = rep(0, dimensions),
                                       Sigma = Sigmaj[names(Sigmaj) == "One"][[1]][1:dimensions, 1:dimensions]),
                       "Two" = mvrnorm(n = SampleSize, mu = rep(0, dimensions),
                                       Sigma = Sigmaj[names(Sigmaj) == "Two"][[1]][1:dimensions, 1:dimensions])),
                  simplify = FALSE)
}, SampleSize = grid$SampleSize, dimensions = grid$dimensions, Sig = grid$Sig,
MoreArgs = list(Sigmaj = Sigmaj, replications = replications), SIMPLIFY = FALSE)

save(mvndata, file = "E:/Ben/Box Sync/Statistics/Unstructured/mvndata.RData")

pushover(message = "mvndata",
         title = "Hey")

NullthreeTests <- ldply(mvndata, function(list){
  ldply(list, function(list){
   data.frame(SampleSize = nrow(list$Zero1), dimension = ncol(list$Zero1),
              Test = c("Chaipitak", "Schott", "Srivastava 2007", "Srivastava 2010",
                       "Srivastava 2014", "Ishii", "Ahmad"),
              Statistic = c(EqualCov:::Chaipitak2013Stat(list(list$Zero1, list$Zero2, list$Zero3)),
                            EqualCov:::Schott2007Stat(list(list$Zero1, list$Zero2, list$Zero3)),
                            EqualCov:::Srivastava2007Stat(list(list$Zero1, list$Zero2, list$Zero3)),
                            EqualCov:::SrivastavaYanagihara2010Stat(list(list$Zero1, list$Zero2, list$Zero3)),
                            EqualCov:::Srivastava2014Stat(list(list$Zero1, list$Zero2, list$Zero3)),
                            EqualCov:::Ishii2016Stat(list(list$Zero1, list$Zero2, list$Zero3)),
                            EqualCov:::Ahmad2017Stat(list(list$Zero1, list$Zero2, list$Zero3))))
  })
})

save(NullthreeTests, file = "E:/Ben/Box Sync/Statistics/Unstructured/NullthreeTests.RData")

pushover(message = "NullthreeTests",
         title = "Hey")

cvsthree <- summarize(group_by(NullthreeTests, SampleSize, dimension, Test),
                      CriticalValue = quantile(Statistic, 0.95))

save(cvsthree, file = "E:/Ben/Box Sync/Statistics/Unstructured/cvsthree.RData")

pushover(message = "cvsthree",
         title = "Hey")

Powervaluesthreetests <- ldply(mvndata, function(list){
  ldply(list, function(list){
    data.frame(SampleSize = nrow(list$Zero1), dimension = ncol(list$Zero1),
               Test = c("Chaipitak", "Schott", "Srivastava 2007", "Srivastava 2010",
                        "Srivastava 2014", "Ishii", "Ahmad"),
               Statistic = c(EqualCov:::Chaipitak2013Stat(list(list$Zero1, list$One, list$Two)),
                             EqualCov:::Schott2007Stat(list(list$Zero1, list$One, list$Two)),
                             EqualCov:::Srivastava2007Stat(list(list$Zero1, list$One, list$Two)),
                             EqualCov:::SrivastavaYanagihara2010Stat(list(list$Zero1, list$One, list$Two)),
                             EqualCov:::Srivastava2014Stat(list(list$Zero1, list$One, list$Two)),
                             EqualCov:::Ishii2016Stat(list(list$Zero1, list$One, list$Two)),
                             EqualCov:::Ahamd2017Stat(list(list$Zero1, list$One, list$Two))))
  })
})

save(Powervaluesthreetests, file = "E:/Ben/Box Sync/Statistics/Unstructured/Powervaluesthreetests.RData")

pushover(message = "Powervaluesthreetests",
         title = "Hey")

powerscoresthreetests <- mutate(full_join(cvsthree, Powervaluesthreetests),
                                Significant = (Statistic > CriticalValue))

save(powerscoresthreetests, file = "E:/Ben/Box Sync/Statistics/Unstructured/powerscoresthreetests.RData")

pushover(message = "powerscoresthreetests",
         title = "Hey")

powerthreetest <- summarise(group_by(powerscoresthreetests,
                                     SampleSize, dimension, Test),
                            Power = mean(Significant))

save(powerthreetest, file = "E:/Ben/Box Sync/Statistics/Unstructured/powerthreetest.RData")

pushover(message = "powerthreetest",
         title = "Hey")

NulltwoTests <- ldply(mvndata, function(list){
  ldply(list, function(list){
    data.frame(SampleSize = nrow(list$Zero1), dimension = ncol(list$Zero1),
               Test = c("Chaipitak", "Schott", "Srivastava 2007", "Srivastava 2010",
                        "Srivastava 2014", "Ishii", "Ahmad"),
               Statistic = c(EqualCov:::Chaipitak2013Stat(list(list$Zero1, list$Zero2)),
                             EqualCov:::Schott2007Stat(list(list$Zero1, list$Zero2)),
                             EqualCov:::Srivastava2007Stat(list(list$Zero1, list$Zero2)),
                             EqualCov:::SrivastavaYanagihara2010Stat(list(list$Zero1, list$Zero2)),
                             EqualCov:::Srivastava2014Stat(list(list$Zero1, list$Zero2)),
                             EqualCov:::Ishii2016Stat(list(list$Zero1, list$Zero2)),
                             EqualCov:::Ahmad2017Stat(list(list$Zero1, list$Zero2))))
  })
})

save(NulltwoTests, file = "E:/Ben/Box Sync/Statistics/Unstructured/NulltwoTests.RData")

pushover(message = "NulltwoTests",
         title = "Hey")

cvstwo <- summarize(group_by(NulltwoTests, SampleSize, dimension, Test),
                      CriticalValue = quantile(Statistic, 0.95))

save(cvstwo, file = "E:/Ben/Box Sync/Statistics/Unstructured/cvstwo.RData")

pushover(message = "cvstwo",
         title = "Hey")

Powervaluestwotests <- ldply(mvndata, function(list){
  ldply(list, function(list){
    data.frame(SampleSize = nrow(list$Zero1), dimension = ncol(list$Zero1),
               Test = c("Chaipitak", "Schott", "Srivastava 2007", "Srivastava 2010",
                        "Srivastava 2014", "Ishii", "Ahmad"),
               Statistic = c(EqualCov:::Chaipitak2013Stat(list(list$Zero1, list$One)),
                             EqualCov:::Schott2007Stat(list(list$Zero1, list$One)),
                             EqualCov:::Srivastava2007Stat(list(list$Zero1, list$One)),
                             EqualCov:::SrivastavaYanagihara2010Stat(list(list$Zero1, list$One)),
                             EqualCov:::Srivastava2014Stat(list(list$Zero1, list$One)),
                             EqualCov:::Ishii2016Stat(list(list$Zero1, list$One)),
                             EqualCov:::Ahmad2017Stat(list(list$Zero1, list$One))))
  })
})

save(Powervaluestwotests, file = "E:/Ben/Box Sync/Statistics/Unstructured/Powervaluestwotests.RData")

pushover(message = "Powervaluestwotests",
         title = "Hey")

powerscorestwotests <- mutate(full_join(cvstwo, Powervaluestwotests),
                              Significant = (Statistic > CriticalValue))

save(powerscorestwotests, file = "E:/Ben/Box Sync/Statistics/Unstructured/powerscorestwotests.RData")

pushover(message = "powerscorestwotests",
         title = "Hey")

powertwotest <- summarise(group_by(powerscorestwotests,
                                   SampleSize, dimension, Test),
                          Power = mean(Significant))

save(powertwotest, file = "E:/Ben/Box Sync/Statistics/Unstructured/powertwotest.RData")

pushover(message = "Your Sim is Finished",
         title = "Hey")
