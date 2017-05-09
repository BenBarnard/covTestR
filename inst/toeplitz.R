library(MASS)
library(plyr)
library(EqualCov)
library(dplyr)
library(pushoverr)
set_pushover_user(user = "ufmfa6vc9s2fc2eh6phop9ej5ebxum")
set_pushover_app(token = "azrd3hwrwgh2gs6igbvb8yy4mftoi7")

dimensions <- c(20, 40, 60, 100, 200)
SampleSize <- c(10, 20, 40, 60)
replications <- 1000
grid <- expand.grid(dimensions = dimensions, SampleSize = SampleSize, Sigmaj = c("Zero", "One", "Two"))

Sigmaj <- list("Zero" = .99 * diag(1, 200) + 0.01 * rep(1, 200) %*% t(rep(1, 200)),
               "One" = .95 * diag(1, 200) + 0.05 * rep(1, 200) %*% t(rep(1, 200)),
               "Two" = .99 * diag(1, 200) + 0.01 * rep(1, 200) %*% t(rep(1, 200)))

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

save(mvndata, file = "E:/Ben/Box Sync/Statistics/toeplitz/mvndata.RData")

pushover(message = "mvndata",
         title = "Hey")

NullthreeTests <- ldply(mvndata, function(list){
  ldply(list, function(list){
    data.frame(SampleSize = nrow(list$Zero1), dimension = ncol(list$Zero1),
               Test = c("Chaipitak", "Schott", "Schott Sample Cov", "Srivastava 2007", "Srivastava 2010",
                        "Srivastava 2014", "Ishii"),
               Statistic = c(Chaipitak2013(list(list$Zero1, list$Zero2, list$Zero3))$statistic,
                             Schott2007(list(list$Zero1, list$Zero2, list$Zero3))$statistic,
                             Schott2007sample(list(list$Zero1, list$Zero2, list$Zero3))$statistic,
                             Srivastava2007(list(list$Zero1, list$Zero2, list$Zero3))$statistic,
                             SrivastavaYanagihara2010(list(list$Zero1, list$Zero2, list$Zero3))$statistic,
                             Srivastava2014(list(list$Zero1, list$Zero2, list$Zero3))$statistic,
                             Ishii2016(list(list$Zero1, list$Zero2, list$Zero3))$statistic))
  })
})

save(NullthreeTests, file = "E:/Ben/Box Sync/Statistics/toeplitz/NullthreeTests.RData")

pushover(message = "NullthreeTests",
         title = "Hey")

cvsthree <- summarize(group_by(NullthreeTests, SampleSize, dimension, Test),
                      CriticalValue = quantile(Statistic, 0.95))

save(cvsthree, file = "E:/Ben/Box Sync/Statistics/toeplitz/cvsthree.RData")

pushover(message = "cvsthree",
         title = "Hey")

Powervaluesthreetests <- ldply(mvndata, function(list){
  ldply(list, function(list){
    data.frame(SampleSize = nrow(list$Zero1), dimension = ncol(list$Zero1),
               Test = c("Chaipitak", "Schott", "Schott Sample Cov", "Srivastava 2007", "Srivastava 2010",
                        "Srivastava 2014", "Ishii"),
               Statistic = c(Chaipitak2013(list(list$Zero1, list$One, list$Two))$statistic,
                             Schott2007(list(list$Zero1, list$One, list$Two))$statistic,
                             Schott2007sample(list(list$Zero1, list$One, list$Two))$statistic,
                             Srivastava2007(list(list$Zero1, list$One, list$Two))$statistic,
                             SrivastavaYanagihara2010(list(list$Zero1, list$One, list$Two))$statistic,
                             Srivastava2014(list(list$Zero1, list$One, list$Two))$statistic,
                             Ishii2016(list(list$Zero1, list$One, list$Two))$statistic))
  })
})

save(Powervaluesthreetests, file = "E:/Ben/Box Sync/Statistics/toeplitz/Powervaluesthreetests.RData")

pushover(message = "Powervaluesthreetests",
         title = "Hey")

powerscoresthreetests <- mutate(full_join(cvsthree, Powervaluesthreetests),
                                Significant = (Statistic > CriticalValue))

save(powerscoresthreetests, file = "E:/Ben/Box Sync/Statistics/toeplitz/powerscoresthreetests.RData")

pushover(message = "powerscoresthreetests",
         title = "Hey")

powerthreetest <- summarise(group_by(powerscoresthreetests,
                                     SampleSize, dimension, Test),
                            Power = mean(Significant))

save(powerthreetest, file = "E:/Ben/Box Sync/Statistics/toeplitz/powerthreetest.RData")

pushover(message = "powerthreetest",
         title = "Hey")

NulltwoTests <- ldply(mvndata, function(list){
  ldply(list, function(list){
    data.frame(SampleSize = nrow(list$Zero1), dimension = ncol(list$Zero1),
               Test = c("Chaipitak", "Schott", "Schott Sample Cov", "Srivastava 2007", "Srivastava 2010",
                        "Srivastava 2014", "Ishii"),
               Statistic = c(Chaipitak2013(list(list$Zero1, list$Zero2))$statistic,
                             Schott2007(list(list$Zero1, list$Zero2))$statistic,
                             Schott2007sample(list(list$Zero1, list$Zero2))$statistic,
                             Srivastava2007(list(list$Zero1, list$Zero2))$statistic,
                             SrivastavaYanagihara2010(list(list$Zero1, list$Zero2))$statistic,
                             Srivastava2014(list(list$Zero1, list$Zero2))$statistic,
                             Ishii2016(list(list$Zero1, list$Zero2))$statistic))
  })
})

save(NulltwoTests, file = "E:/Ben/Box Sync/Statistics/toeplitz/NulltwoTests.RData")

pushover(message = "NulltwoTests",
         title = "Hey")

cvstwo <- summarize(group_by(NulltwoTests, SampleSize, dimension, Test),
                    CriticalValue = quantile(Statistic, 0.95))

save(cvstwo, file = "E:/Ben/Box Sync/Statistics/toeplitz/cvstwo.RData")

pushover(message = "cvstwo",
         title = "Hey")

Powervaluestwotests <- ldply(mvndata, function(list){
  ldply(list, function(list){
    data.frame(SampleSize = nrow(list$Zero1), dimension = ncol(list$Zero1),
               Test = c("Chaipitak", "Schott", "Schott Sample Cov", "Srivastava 2007", "Srivastava 2010",
                        "Srivastava 2014", "Ishii"),
               Statistic = c(Chaipitak2013(list(list$Zero1, list$One))$statistic,
                             Schott2007(list(list$Zero1, list$One))$statistic,
                             Schott2007sample(list(list$Zero1, list$One))$statistic,
                             Srivastava2007(list(list$Zero1, list$One))$statistic,
                             SrivastavaYanagihara2010(list(list$Zero1, list$One))$statistic,
                             Srivastava2014(list(list$Zero1, list$One))$statistic,
                             Ishii2016(list(list$Zero1, list$One))$statistic))
  })
})

save(Powervaluestwotests, file = "E:/Ben/Box Sync/Statistics/toeplitz/Powervaluestwotests.RData")

pushover(message = "Powervaluestwotests",
         title = "Hey")

powerscorestwotests <- mutate(full_join(cvstwo, Powervaluestwotests),
                              Significant = (Statistic > CriticalValue))

save(powerscorestwotests, file = "E:/Ben/Box Sync/Statistics/toeplitz/powerscorestwotests.RData")

pushover(message = "powerscorestwotests",
         title = "Hey")

powertwotest <- summarise(group_by(powerscorestwotests,
                                   SampleSize, dimension, Test),
                          Power = mean(Significant))

save(powertwotest, file = "E:/Ben/Box Sync/Statistics/toeplitz/powertwotest.RData")

pushover(message = "Your Sim is Finished",
         title = "Hey")
