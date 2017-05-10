library(MASS)
library(plyr)
library(EqualCov)
library(dplyr)
library(pushoverr)
set_pushover_user(user = "ufmfa6vc9s2fc2eh6phop9ej5ebxum")
set_pushover_app(token = "azrd3hwrwgh2gs6igbvb8yy4mftoi7")

dimensions <- 200
SampleSize <- 10
replications <- 1000

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

save(Sigmaj, file = "E:/Ben/Box Sync/Statistics/Srivastava2010SimDimReduce/Sigmaj.RData")

mvndata <- replicate(replications,
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

save(mvndata, file = "E:/Ben/Box Sync/Statistics/Srivastava2010SimDimReduce/mvndata.RData")

pushover(message = "mvndata",
         title = "Hey")

NullthreeTests <- ldply(mvndata, function(ls, SampleSize, dimensions){
  redMat <- songEquality(ls[1:3])$u
  ldply(c(1:SampleSize, dimensions), function(reduction, ls, redMat){
    lt <- ls
    SampleSize <- nrow(lt$Zero1)
    originaldimension <- ncol(lt$Zero1)
    projection <- t(redMat[,1:reduction])

    if(reduction == originaldimension){
      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                       reduction = reduction,
                       Test = c("Chaipitak", "Schott Sample Cov", "Schott", "Srivastava 2007", "Srivastava 2010",
                                "Srivastava 2014"),
                       Statistic = c(Chaipitak2013(list(lt$Zero1, lt$Zero2, lt$Zero3))$statistic,
                                     Schott2007sample(list(lt$Zero1, lt$Zero2, lt$Zero3))$statistic,
                                     Schott2007(list(lt$Zero1, lt$Zero2, lt$Zero3))$statistic,
                                     Srivastava2007(list(lt$Zero1, lt$Zero2, lt$Zero3))$statistic,
                                     SrivastavaYanagihara2010(list(lt$Zero1, lt$Zero2, lt$Zero3))$statistic,
                                     Srivastava2014(list(lt$Zero1, lt$Zero2, lt$Zero3))$statistic))

    }

    lt <- lapply(lt, function(x){
      t(projection %*% t(x))
    })

    if(reduction < SampleSize){
      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                 reduction = reduction,
                 Test = c("Schott 2001", "Box", "Chaipitak", "Schott Sample Cov", "Schott", "Srivastava 2007", "Srivastava 2010",
                          "Srivastava 2014"),
                 Statistic = c(Schott2001(list(lt$Zero1, lt$Zero2, lt$Zero3))$statistic,
                               BoxesM(list(lt$Zero1, lt$Zero2, lt$Zero3))$statistic,
                               Chaipitak2013(list(lt$Zero1, lt$Zero2, lt$Zero3))$statistic,
                               Schott2007sample(list(lt$Zero1, lt$Zero2, lt$Zero3))$statistic,
                               Schott2007(list(lt$Zero1, lt$Zero2, lt$Zero3))$statistic,
                               Srivastava2007(list(lt$Zero1, lt$Zero2, lt$Zero3))$statistic,
                               SrivastavaYanagihara2010(list(lt$Zero1, lt$Zero2, lt$Zero3))$statistic,
                               Srivastava2014(list(lt$Zero1, lt$Zero2, lt$Zero3))$statistic))

    }

    if(reduction == SampleSize){
      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                       reduction = reduction,
                       Test = c("Chaipitak", "Schott Sample Cov", "Schott", "Srivastava 2007", "Srivastava 2010",
                                "Srivastava 2014"),
                       Statistic = c(Chaipitak2013(list(lt$Zero1, lt$Zero2, lt$Zero3))$statistic,
                                     Schott2007sample(list(lt$Zero1, lt$Zero2, lt$Zero3))$statistic,
                                     Schott2007(list(lt$Zero1, lt$Zero2, lt$Zero3))$statistic,
                                     Srivastava2007(list(lt$Zero1, lt$Zero2, lt$Zero3))$statistic,
                                     SrivastavaYanagihara2010(list(lt$Zero1, lt$Zero2, lt$Zero3))$statistic,
                                     Srivastava2014(list(lt$Zero1, lt$Zero2, lt$Zero3))$statistic))

    }

    df

  }, ls = ls, redMat = redMat)
}, SampleSize = SampleSize, dimensions = dimensions)

save(NullthreeTests, file = "E:/Ben/Box Sync/Statistics/Srivastava2010SimDimReduce/NullthreeTests.RData")

pushover(message = "NullthreeTests",
         title = "Hey")

cvsthree <- summarize(group_by(NullthreeTests, SampleSize, originaldimension, reduction, Test),
                      CriticalValue = quantile(Statistic, 0.95))

save(cvsthree, file = "E:/Ben/Box Sync/Statistics/Srivastava2010SimDimReduce/cvsthree.RData")

pushover(message = "cvsthree",
         title = "Hey")

Powervaluesthreetests <- ldply(mvndata, function(ls, SampleSize, dimensions){
  redMat <- songEquality(ls[1:3])$u
  ldply(c(1:SampleSize, dimensions), function(reduction, ls, redMat){
    lt <- ls
    SampleSize <- nrow(lt$Zero1)
    originaldimension <- ncol(lt$Zero1)
    projection <- t(redMat[,1:reduction])

    if(reduction == originaldimension){
      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                       reduction = reduction,
                       Test = c("Chaipitak", "Schott Sample Cov", "Schott", "Srivastava 2007", "Srivastava 2010",
                                "Srivastava 2014"),
                       Statistic = c(Chaipitak2013(list(lt$Zero1, lt$One, lt$Two))$statistic,
                                     Schott2007sample(list(lt$Zero1, lt$One, lt$Two))$statistic,
                                     Schott2007(list(lt$Zero1, lt$One, lt$Two))$statistic,
                                     Srivastava2007(list(lt$Zero1, lt$One, lt$Two))$statistic,
                                     SrivastavaYanagihara2010(list(lt$Zero1, lt$One, lt$Two))$statistic,
                                     Srivastava2014(list(lt$Zero1, lt$One, lt$Two))$statistic))

    }

    lt <- lapply(lt, function(x){
      t(projection %*% t(x))
    })

    if(reduction < SampleSize){
      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                       reduction = reduction,
                       Test = c("Schott 2001", "Box", "Chaipitak", "Schott Sample Cov", "Schott", "Srivastava 2007", "Srivastava 2010",
                                "Srivastava 2014"),
                       Statistic = c(Schott2001(list(lt$Zero1, lt$One, lt$Two))$statistic,
                                     BoxesM(list(lt$Zero1, lt$One, lt$Two))$statistic,
                                     Chaipitak2013(list(lt$Zero1, lt$One, lt$Two))$statistic,
                                     Schott2007sample(list(lt$Zero1, lt$One, lt$Two))$statistic,
                                     Schott2007(list(lt$Zero1, lt$One, lt$Two))$statistic,
                                     Srivastava2007(list(lt$Zero1, lt$One, lt$Two))$statistic,
                                     SrivastavaYanagihara2010(list(lt$Zero1, lt$One, lt$Two))$statistic,
                                     Srivastava2014(list(lt$Zero1, lt$One, lt$Two))$statistic))

    }

    if(reduction == SampleSize){
      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                       reduction = reduction,
                       Test = c("Chaipitak", "Schott Sample Cov", "Schott", "Srivastava 2007", "Srivastava 2010",
                                "Srivastava 2014"),
                       Statistic = c(Chaipitak2013(list(lt$Zero1, lt$One, lt$Two))$statistic,
                                     Schott2007sample(list(lt$Zero1, lt$One, lt$Two))$statistic,
                                     Schott2007(list(lt$Zero1, lt$One, lt$Two))$statistic,
                                     Srivastava2007(list(lt$Zero1, lt$One, lt$Two))$statistic,
                                     SrivastavaYanagihara2010(list(lt$Zero1, lt$One, lt$Two))$statistic,
                                     Srivastava2014(list(lt$Zero1, lt$One, lt$Two))$statistic))

    }

    df

  }, ls = ls, redMat = redMat)
}, SampleSize = SampleSize, dimensions = dimensions)

save(Powervaluesthreetests, file = "E:/Ben/Box Sync/Statistics/Srivastava2010SimDimReduce/Powervaluesthreetests.RData")

pushover(message = "Powervaluesthreetests",
         title = "Hey")

powerscoresthreetests <- mutate(full_join(cvsthree, Powervaluesthreetests),
                                Significant = (Statistic > CriticalValue))

save(powerscoresthreetests, file = "E:/Ben/Box Sync/Statistics/Srivastava2010SimDimReduce/powerscoresthreetests.RData")

pushover(message = "powerscoresthreetests",
         title = "Hey")

powerthreetest <- summarise(group_by(powerscoresthreetests,
                                     SampleSize, originaldimension, reduction, Test),
                            Power = mean(Significant))

save(powerthreetest, file = "E:/Ben/Box Sync/Statistics/Srivastava2010SimDimReduce/powerthreetest.RData")

pushover(message = "powerthreetest",
         title = "Hey")

NulltwoTests <- ldply(mvndata, function(ls, SampleSize, dimensions){
  redMat <- songEquality(ls[1:3])$u
  ldply(c(1:SampleSize, dimensions), function(reduction, ls, redMat){
    lt <- ls
    SampleSize <- nrow(lt$Zero1)
    originaldimension <- ncol(lt$Zero1)
    projection <- t(redMat[,1:reduction])

    if(reduction == originaldimension){
      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                       reduction = reduction,
                       Test = c("Chaipitak", "Schott Sample Cov", "Schott", "Srivastava 2007", "Srivastava 2010",
                                "Srivastava 2014"),
                       Statistic = c(Chaipitak2013(list(lt$Zero1, lt$Zero2))$statistic,
                                     Schott2007sample(list(lt$Zero1, lt$Zero2))$statistic,
                                     Schott2007(list(lt$Zero1, lt$Zero2))$statistic,
                                     Srivastava2007(list(lt$Zero1, lt$Zero2))$statistic,
                                     SrivastavaYanagihara2010(list(lt$Zero1, lt$Zero2))$statistic,
                                     Srivastava2014(list(lt$Zero1, lt$Zero2))$statistic))

    }

    lt <- lapply(lt, function(x){
      t(projection %*% t(x))
    })

    if(reduction < SampleSize){
      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                       reduction = reduction,
                       Test = c("Schott 2001", "Box", "Chaipitak", "Schott Sample Cov", "Schott", "Srivastava 2007", "Srivastava 2010",
                                "Srivastava 2014"),
                       Statistic = c(Schott2001(list(lt$Zero1, lt$Zero2))$statistic,
                                     BoxesM(list(lt$Zero1, lt$Zero2))$statistic,
                                     Chaipitak2013(list(lt$Zero1, lt$Zero2))$statistic,
                                     Schott2007sample(list(lt$Zero1, lt$Zero2))$statistic,
                                     Schott2007(list(lt$Zero1, lt$Zero2))$statistic,
                                     Srivastava2007(list(lt$Zero1, lt$Zero2))$statistic,
                                     SrivastavaYanagihara2010(list(lt$Zero1, lt$Zero2))$statistic,
                                     Srivastava2014(list(lt$Zero1, lt$Zero2))$statistic))

    }

    if(reduction == SampleSize){
      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                       reduction = reduction,
                       Test = c("Chaipitak", "Schott Sample Cov", "Schott", "Srivastava 2007", "Srivastava 2010",
                                "Srivastava 2014"),
                       Statistic = c(Chaipitak2013(list(lt$Zero1, lt$Zero2))$statistic,
                                     Schott2007sample(list(lt$Zero1, lt$Zero2))$statistic,
                                     Schott2007(list(lt$Zero1, lt$Zero2))$statistic,
                                     Srivastava2007(list(lt$Zero1, lt$Zero2))$statistic,
                                     SrivastavaYanagihara2010(list(lt$Zero1, lt$Zero2))$statistic,
                                     Srivastava2014(list(lt$Zero1, lt$Zero2))$statistic))

    }

    df

  }, ls = ls, redMat = redMat)
}, SampleSize = SampleSize, dimensions = dimensions)

save(NulltwoTests, file = "E:/Ben/Box Sync/Statistics/Srivastava2010SimDimReduce/NulltwoTests.RData")

pushover(message = "NulltwoTests",
         title = "Hey")

cvstwo <- summarize(group_by(NulltwoTests, SampleSize, originaldimension, reduction, Test),
                      CriticalValue = quantile(Statistic, 0.95))

save(cvstwo, file = "E:/Ben/Box Sync/Statistics/Srivastava2010SimDimReduce/cvstwo.RData")

pushover(message = "cvstwo",
         title = "Hey")

Powervaluestwotests <- ldply(mvndata, function(ls, SampleSize, dimensions){
  redMat <- songEquality(ls[1:3])$u
  ldply(c(1:SampleSize, dimensions), function(reduction, ls, redMat){
    lt <- ls
    SampleSize <- nrow(lt$Zero1)
    originaldimension <- ncol(lt$Zero1)
    projection <- t(redMat[,1:reduction])

    if(reduction == originaldimension){
      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                       reduction = reduction,
                       Test = c("Chaipitak", "Schott Sample Cov", "Schott", "Srivastava 2007", "Srivastava 2010",
                                "Srivastava 2014"),
                       Statistic = c(Chaipitak2013(list(lt$Zero1, lt$One))$statistic,
                                     Schott2007sample(list(lt$Zero1, lt$One))$statistic,
                                     Schott2007(list(lt$Zero1, lt$One))$statistic,
                                     Srivastava2007(list(lt$Zero1, lt$One))$statistic,
                                     SrivastavaYanagihara2010(list(lt$Zero1, lt$One))$statistic,
                                     Srivastava2014(list(lt$Zero1, lt$One))$statistic))

    }

    lt <- lapply(lt, function(x){
      t(projection %*% t(x))
    })

    if(reduction < SampleSize){
      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                       reduction = reduction,
                       Test = c("Schott 2001", "Box", "Chaipitak", "Schott Sample Cov", "Schott", "Srivastava 2007", "Srivastava 2010",
                                "Srivastava 2014"),
                       Statistic = c(Schott2001(list(lt$Zero1, lt$One))$statistic,
                                     BoxesM(list(lt$Zero1, lt$One))$statistic,
                                     Chaipitak2013(list(lt$Zero1, lt$One))$statistic,
                                     Schott2007sample(list(lt$Zero1, lt$One))$statistic,
                                     Schott2007(list(lt$Zero1, lt$One))$statistic,
                                     Srivastava2007(list(lt$Zero1, lt$One))$statistic,
                                     SrivastavaYanagihara2010(list(lt$Zero1, lt$One))$statistic,
                                     Srivastava2014(list(lt$Zero1, lt$One))$statistic))

    }

    if(reduction == SampleSize){
      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                       reduction = reduction,
                       Test = c("Chaipitak", "Schott Sample Cov", "Schott", "Srivastava 2007", "Srivastava 2010",
                                "Srivastava 2014"),
                       Statistic = c(Chaipitak2013(list(lt$Zero1, lt$One))$statistic,
                                     Schott2007sample(list(lt$Zero1, lt$One))$statistic,
                                     Schott2007(list(lt$Zero1, lt$One))$statistic,
                                     Srivastava2007(list(lt$Zero1, lt$One))$statistic,
                                     SrivastavaYanagihara2010(list(lt$Zero1, lt$One))$statistic,
                                     Srivastava2014(list(lt$Zero1, lt$One))$statistic))

    }

    df

  }, ls = ls, redMat = redMat)
}, SampleSize = SampleSize, dimensions = dimensions)

save(Powervaluestwotests, file = "E:/Ben/Box Sync/Statistics/Srivastava2010SimDimReduce/Powervaluestwotests.RData")

pushover(message = "Powervaluestwotests",
         title = "Hey")

powerscorestwotests <- mutate(full_join(cvstwo, Powervaluestwotests),
                              Significant = (Statistic > CriticalValue))

save(powerscorestwotests, file = "E:/Ben/Box Sync/Statistics/Srivastava2010SimDimReduce/powerscorestwotests.RData")

pushover(message = "powerscorestwotests",
         title = "Hey")

powertwotest <- summarise(group_by(powerscorestwotests,
                                   SampleSize, originaldimension, reduction, Test),
                          Power = mean(Significant))

save(powertwotest, file = "E:/Ben/Box Sync/Statistics/Srivastava2010SimDimReduce/powertwotest.RData")

pushover(message = "Your Sim is Finished",
         title = "Hey")
