library(MASS)
library(plyr)
library(EqualCov)
library(dplyr)
library(pushoverr)
set_pushover_user(user = "ufmfa6vc9s2fc2eh6phop9ej5ebxum")
set_pushover_app(token = "azrd3hwrwgh2gs6igbvb8yy4mftoi7")

dimensions <- 20
SampleSize <- 10

load(file = "E:/Ben/Box Sync/Statistics/Autoregressive/mvndata.RData")

NullTests <- ldply(mvndata, function(ls, SampleSize, dimensions){
  ls <- lapply(ls, function(x){
    browser()
    x[1:SampleSize, 1:dimensions]
  })[1:3]
  covs <- lapply(ls, cov)
  diffs <- lapply(covs, function(x){x - covs[[1]]})[-1]
  redMat2 <- svd(Reduce(cbind, diffs[-2]))$u
  redMat3 <- svd(Reduce(cbind, diffs))$u

  ldply(c(1:SampleSize, dimensions), function(reduction, ls, redMat2, redMat3, SampleSize, dimensions){
    lt <- ls
    projection2 <- t(redMat2[,1:reduction])
    projection3 <- t(redMat3[,1:reduction])
    SampleSize <- SampleSize
    originaldimension <- dimensions

    if(reduction == originaldimension){
      lsmat <- list(lt$Zero1, lt$Zero2, lt$Zero3)
      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                       reduction = reduction, Pops = c(rep("Two", 6), rep("Three", 6)),
                       Test = c("Chaipitak", "Ahmad", "Schott", "Srivastava 2007", "Srivastava 2010",
                                "Srivastava 2014", "Chaipitak", "Ahmad", "Schott", "Srivastava 2007",
                                "Srivastava 2010", "Srivastava 2014"),
                       Statistic = c(EqualCov:::Chaipitak2013Stat(lsmat),
                                     EqualCov:::Ahmad2017Stat(lsmat),
                                     EqualCov:::Schott2007Stat(lsmat),
                                     EqualCov:::Srivastava2007Stat(lsmat),
                                     EqualCov:::SrivastavaYanagihara2010Stat(lsmat),
                                     EqualCov:::Srivastava2014Stat(lsmat),
                                     EqualCov:::Chaipitak2013Stat(lsmat[-3]),
                                     EqualCov:::Ahmad2017Stat(lsmat[-3]),
                                     EqualCov:::Schott2007Stat(lsmat[-3]),
                                     EqualCov:::Srivastava2007Stat(lsmat[-3]),
                                     EqualCov:::SrivastavaYanagihara2010Stat(lsmat[-3]),
                                     EqualCov:::Srivastava2014Stat(lsmat[-3])))

    }

    lt2 <- lapply(lt, function(x){
      t(projection2 %*% t(x))
    })
    lt3 <- lapply(lt, function(x){
      t(projection3 %*% t(x))
    })
    lsmat3 <- list(lt3$Zero1, lt3$Zero2, lt3$Zero3)
    lsmat2 <- list(lt2$Zero1, lt2$Zero2)

    if(reduction < SampleSize){
      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                       reduction = reduction, Pops = c(rep("Two", 8), rep("Three", 8)),
                       Test = c("Schott 2001", "Box","Chaipitak", "Ahmad", "Schott", "Srivastava 2007",
                                "Srivastava 2010", "Srivastava 2014", "Schott 2001", "Box","Chaipitak",
                                "Ahmad", "Schott", "Srivastava 2007", "Srivastava 2010", "Srivastava 2014"),
                       Statistic = c(EqualCov:::Schott2001Stat(lsmat3),
                                     EqualCov:::BoxesMStat(lsmat3),
                                     EqualCov:::Chaipitak2013Stat(lsmat3),
                                     EqualCov:::Ahmad2017Stat(lsmat3),
                                     EqualCov:::Schott2007Stat(lsmat3),
                                     EqualCov:::Srivastava2007Stat(lsmat3),
                                     EqualCov:::SrivastavaYanagihara2010Stat(lsmat3),
                                     EqualCov:::Srivastava2014Stat(lsmat3),
                                     EqualCov:::Schott2001Stat(lsmat2),
                                     EqualCov:::BoxesMStat(lsmat2),
                                     EqualCov:::Chaipitak2013Stat(lsmat2),
                                     EqualCov:::Ahmad2017Stat(lsmat2),
                                     EqualCov:::Schott2007Stat(lsmat2),
                                     EqualCov:::Srivastava2007Stat(lsmat2),
                                     EqualCov:::SrivastavaYanagihara2010Stat(lsmat2),
                                     EqualCov:::Srivastava2014Stat(lsmat2)))

    }

    if(reduction == SampleSize){
      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                       reduction = reduction, Pops = c(rep("Two", 6), rep("Three", 6)),
                       Test = c("Chaipitak", "Ahmad", "Schott", "Srivastava 2007", "Srivastava 2010",
                                "Srivastava 2014", "Chaipitak", "Ahmad", "Schott", "Srivastava 2007",
                                "Srivastava 2010", "Srivastava 2014"),
                       Statistic = c(EqualCov:::Chaipitak2013Stat(lsmat3),
                                     EqualCov:::Ahmad2017Stat(lsmat3),
                                     EqualCov:::Schott2007Stat(lsmat3),
                                     EqualCov:::Srivastava2007Stat(lsmat3),
                                     EqualCov:::SrivastavaYanagihara2010Stat(lsmat3),
                                     EqualCov:::Srivastava2014Stat(lsmat3),
                                     EqualCov:::Chaipitak2013Stat(lsmat2),
                                     EqualCov:::Ahmad2017Stat(lsmat2),
                                     EqualCov:::Schott2007Stat(lsmat2),
                                     EqualCov:::Srivastava2007Stat(lsmat2),
                                     EqualCov:::SrivastavaYanagihara2010Stat(lsmat2),
                                     EqualCov:::Srivastava2014Stat(lsmat2)))

    }

    df

  }, ls = ls, redMat2 = redMat2, redMat3 = redMat3, SampleSize = SampleSize, dimensions = dimensions)
}, SampleSize = SampleSize, dimensions = dimensions)

save(NullTests, file = "E:/Ben/Box Sync/Statistics/Autoregressive/DimReduce/NullTests.RData")

pushover(message = "NullTests",
         title = "Hey")

cvsthree <- summarize(group_by(NullTests, SampleSize, originaldimension, reduction, Test),
                      CriticalValue = quantile(Statistic, 0.95))

save(cvs, file = "E:/Ben/Box Sync/Statistics/Autoregressive/DimReduce/cvs.RData")

pushover(message = "cvs",
         title = "Hey")

Powervaluestests <- ldply(mvndata, function(ls, SampleSize, dimensions){
  ls <- lapply(ls, function(x){
    x[1:SampleSize, 1:dimensions]
  })[c(1, 4, 5)]
  covs <- lapply(ls, cov)
  diffs <- lapply(covs, function(x){x - covs[[1]]})[-1]
  redMat2 <- svd(Reduce(cbind, diffs[-2]))$u
  redMat3 <- svd(Reduce(cbind, diffs))$u

  ldply(c(1:SampleSize, dimensions), function(reduction, ls, redMat2, redMat3, SampleSize, dimensions){
    lt <- ls
    projection2 <- t(redMat2[,1:reduction])
    projection3 <- t(redMat3[,1:reduction])
    SampleSize <- SampleSize
    originaldimension <- dimensions

    if(reduction == originaldimension){
      lsmat <- list(lt$Zero1, lt$One, lt$Two)
      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                       reduction = reduction, Pops = c(rep("Two", 6), rep("Three", 6)),
                       Test = c("Chaipitak", "Ahmad", "Schott", "Srivastava 2007", "Srivastava 2010",
                                "Srivastava 2014", "Chaipitak", "Ahmad", "Schott", "Srivastava 2007",
                                "Srivastava 2010", "Srivastava 2014"),
                       Statistic = c(EqualCov:::Chaipitak2013Stat(lsmat),
                                     EqualCov:::Ahmad2017Stat(lsmat),
                                     EqualCov:::Schott2007Stat(lsmat),
                                     EqualCov:::Srivastava2007Stat(lsmat),
                                     EqualCov:::SrivastavaYanagihara2010Stat(lsmat),
                                     EqualCov:::Srivastava2014Stat(lsmat),
                                     EqualCov:::Chaipitak2013Stat(lsmat[-3]),
                                     EqualCov:::Ahmad2017Stat(lsmat[-3]),
                                     EqualCov:::Schott2007Stat(lsmat[-3]),
                                     EqualCov:::Srivastava2007Stat(lsmat[-3]),
                                     EqualCov:::SrivastavaYanagihara2010Stat(lsmat[-3]),
                                     EqualCov:::Srivastava2014Stat(lsmat[-3])))

    }

    lt2 <- lapply(lt, function(x){
      t(projection2 %*% t(x))
    })
    lt3 <- lapply(lt, function(x){
      t(projection3 %*% t(x))
    })
    lsmat3 <- list(lt3$Zero1, lt3$One, lt3$Two)
    lsmat2 <- list(lt2$Zero1, lt2$One)

    if(reduction < SampleSize){

      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                       reduction = reduction, Pops = c(rep("Two", 8), rep("Three", 8)),
                       Test = c("Schott 2001", "Box","Chaipitak", "Ahmad", "Schott", "Srivastava 2007",
                                "Srivastava 2010", "Srivastava 2014", "Schott 2001", "Box","Chaipitak",
                                "Ahmad", "Schott", "Srivastava 2007", "Srivastava 2010", "Srivastava 2014"),
                       Statistic = c(EqualCov:::Schott2001Stat(lsmat3),
                                     EqualCov:::BoxesMStat(lsmat3),
                                     EqualCov:::Chaipitak2013Stat(lsmat3),
                                     EqualCov:::Ahmad2017Stat(lsmat3),
                                     EqualCov:::Schott2007Stat(lsmat3),
                                     EqualCov:::Srivastava2007Stat(lsmat3),
                                     EqualCov:::SrivastavaYanagihara2010Stat(lsmat3),
                                     EqualCov:::Srivastava2014Stat(lsmat3),
                                     EqualCov:::Schott2001Stat(lsmat2),
                                     EqualCov:::BoxesMStat(lsmat2),
                                     EqualCov:::Chaipitak2013Stat(lsmat2),
                                     EqualCov:::Ahmad2017Stat(lsmat2),
                                     EqualCov:::Schott2007Stat(lsmat2),
                                     EqualCov:::Srivastava2007Stat(lsmat2),
                                     EqualCov:::SrivastavaYanagihara2010Stat(lsmat2),
                                     EqualCov:::Srivastava2014Stat(lsmat2)))

    }

    if(reduction == SampleSize){

      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                       reduction = reduction, Pops = c(rep("Two", 6), rep("Three", 6)),
                       Test = c("Chaipitak", "Ahmad", "Schott", "Srivastava 2007", "Srivastava 2010",
                                "Srivastava 2014", "Chaipitak", "Ahmad", "Schott", "Srivastava 2007",
                                "Srivastava 2010", "Srivastava 2014"),
                       Statistic = c(EqualCov:::Chaipitak2013Stat(lsmat3),
                                     EqualCov:::Ahmad2017Stat(lsmat3),
                                     EqualCov:::Schott2007Stat(lsmat3),
                                     EqualCov:::Srivastava2007Stat(lsmat3),
                                     EqualCov:::SrivastavaYanagihara2010Stat(lsmat3),
                                     EqualCov:::Srivastava2014Stat(lsmat3),
                                     EqualCov:::Chaipitak2013Stat(lsmat2),
                                     EqualCov:::Ahmad2017Stat(lsmat2),
                                     EqualCov:::Schott2007Stat(lsmat2),
                                     EqualCov:::Srivastava2007Stat(lsmat2),
                                     EqualCov:::SrivastavaYanagihara2010Stat(lsmat2),
                                     EqualCov:::Srivastava2014Stat(lsmat2)))

    }

    df

  }, ls = ls, redMat2 = redMat2, redMat3 = redMat3, SampleSize = SampleSize, dimensions = dimensions)
}, SampleSize = SampleSize, dimensions = dimensions)

save(Powervaluestests, file = "E:/Ben/Box Sync/Statistics/Autoregressive/DimReduce/Powervaluestests.RData")

pushover(message = "Powervaluestests",
         title = "Hey")

powerscorestests <- mutate(full_join(cvs, Powervaluestests),
                                Significant = (Statistic > CriticalValue))

save(powerscorestests, file = "E:/Ben/Box Sync/Statistics/Autoregressive/DimReduce/powerscorestests.RData")

pushover(message = "powerscorestests",
         title = "Hey")

power <- summarise(group_by(powerscorestests,
                                     SampleSize, originaldimension, reduction, Test),
                            Power = mean(Significant))

save(powerthreetest, file = "E:/Ben/Box Sync/Statistics/Autoregressive/DimReduce/power.RData")

pushover(message = "power",
         title = "Hey")
