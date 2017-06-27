source("E:/Ben/EqualCov/inst/reductionFunctions.R")
library(MASS)
library(plyr)
library(EqualCov)
library(dplyr)
library(pushoverr)
set_pushover_user(user = "ufmfa6vc9s2fc2eh6phop9ej5ebxum")
set_pushover_app(token = "azrd3hwrwgh2gs6igbvb8yy4mftoi7")

dimensions <- 20
SampleSize <- 10
gridcomb <- filter(expand.grid(Samples = SampleSize,
                               dims = dimensions),
                   dims > Samples)
reductionMethod <- sdiff

Structure <- "Sri2014"

if(!(dir.exists(paste0("E:/Ben/Box Sync/Statistics/Dissertation/Sims/", Structure, "/DimReduce/", dimensions, " ", SampleSize, "/")))){
  dir.create(paste0("E:/Ben/Box Sync/Statistics/Dissertation/Sims/", Structure, "/DimReduce/", dimensions, " ", SampleSize, "/"))
}


load(file = paste0("E:/Ben/Box Sync/Statistics/Dissertation/mvnData/", Structure, "/mvndata.RData"))

statistics <- ldply(mvndata, function(ls, SampleSize, dimensions, reductionMethod){

    nullList <- list(ls$Zero1[1:SampleSize, 1:dimensions],
                     ls$Zero2[1:SampleSize, 1:dimensions],
                     ls$Zero3[1:SampleSize, 1:dimensions])

    powerList <- list(ls$Zero1[1:SampleSize, 1:dimensions],
                      ls$One[1:SampleSize, 1:dimensions],
                      ls$Two[1:SampleSize, 1:dimensions])

    nullMat2 <- do.call(reductionMethod, list(nullList[-3]))
    powerMat2 <- do.call(reductionMethod, list(nullList[-3]))
    nullMat3 <- do.call(reductionMethod, list(powerList))
    powerMat3 <- do.call(reductionMethod, list(powerList))

  ldply(c(1:SampleSize, dimensions), function(reduction, nullList, powerList, nullredMat2, nullredMat3,
                                              powerredMat2, powerredMat3, SampleSize, dimensions){


    nullprojection2 <- t(nullredMat2[,1:reduction])
    nullprojection3 <- t(nullredMat3[,1:reduction])
    powerprojection2 <- t(powerredMat2[,1:reduction])
    powerprojection3 <- t(powerredMat3[,1:reduction])


    SampleSize <- SampleSize
    originaldimension <- dimensions

    nullList2 <- nullList[-3]
    powerList2 <- powerList[-3]
    nullList3 <- nullList
    powerList3 <- powerList

    if(reduction < originaldimension){
      nullList2 <- lapply(nullList2, function(x){
        t(nullprojection2 %*% t(x))
      })

      nullList3 <- lapply(nullList3, function(x){
        t(nullprojection3 %*% t(x))
      })

      powerList2 <- lapply(powerList2, function(x){
        t(powerprojection2 %*% t(x))
      })

      powerList3 <- lapply(powerList3, function(x){
        t(powerprojection3 %*% t(x))
      })

    }
      dfTest <- data_frame(SampleSize = SampleSize, dimension = originaldimension, reduction = reduction,
                           Test = rep(c("Chaipitak", "Chaipitak pool", "Schott", "Schott pool",
                                        "Srivastava 2007", "Srivastava 2010", "Srivastava 2014",
                                        "Srivastava 2014 pool", "Ishii", "Ahmad"), 2),
                           Pops = c(rep("Three", 10), rep("Two", 10)),
                           `Null Statistic` = c(EqualCov:::Chaipitak2013Stat(nullList3),
                                                EqualCov:::Chaipitak2013poolStat(nullList3),
                                                EqualCov:::Schott2007Stat(nullList3),
                                                EqualCov:::Schott2007pooledStat(nullList3),
                                                EqualCov:::Srivastava2007Stat(nullList3),
                                                EqualCov:::SrivastavaYanagihara2010Stat(nullList3),
                                                EqualCov:::Srivastava2014Stat(nullList3),
                                                EqualCov:::Srivastava2014poolStat(nullList3),
                                                EqualCov:::Ishii2016Stat(nullList3),
                                                abs(EqualCov:::Ahmad2017Stat(nullList3)),
                                                EqualCov:::Chaipitak2013Stat(nullList2),
                                                EqualCov:::Chaipitak2013poolStat(nullList2),
                                                EqualCov:::Schott2007Stat(nullList2),
                                                EqualCov:::Schott2007pooledStat(nullList2),
                                                EqualCov:::Srivastava2007Stat(nullList2),
                                                EqualCov:::SrivastavaYanagihara2010Stat(nullList2),
                                                EqualCov:::Srivastava2014Stat(nullList2),
                                                EqualCov:::Srivastava2014poolStat(nullList2),
                                                EqualCov:::Ishii2016Stat(nullList2),
                                                abs(EqualCov:::Ahmad2017Stat(nullList2))),
                           `Power Statistic` = c(EqualCov:::Chaipitak2013Stat(powerList3),
                                                 EqualCov:::Chaipitak2013poolStat(powerList3),
                                                 EqualCov:::Schott2007Stat(powerList3),
                                                 EqualCov:::Schott2007pooledStat(powerList3),
                                                 EqualCov:::Srivastava2007Stat(powerList3),
                                                 EqualCov:::SrivastavaYanagihara2010Stat(powerList3),
                                                 EqualCov:::Srivastava2014Stat(powerList3),
                                                 EqualCov:::Srivastava2014poolStat(powerList3),
                                                 EqualCov:::Ishii2016Stat(powerList3),
                                                 abs(EqualCov:::Ahmad2017Stat(powerList3)),
                                                 EqualCov:::Chaipitak2013Stat(powerList2),
                                                 EqualCov:::Chaipitak2013poolStat(powerList2),
                                                 EqualCov:::Schott2007Stat(powerList2),
                                                 EqualCov:::Schott2007pooledStat(powerList2),
                                                 EqualCov:::Srivastava2007Stat(powerList2),
                                                 EqualCov:::SrivastavaYanagihara2010Stat(powerList2),
                                                 EqualCov:::Srivastava2014Stat(powerList2),
                                                 EqualCov:::Srivastava2014poolStat(powerList2),
                                                 EqualCov:::Ishii2016Stat(powerList2),
                                                 abs(EqualCov:::Ahmad2017Stat(powerList2))))

    dfTest

  }, nullList = nullList, powerList = powerList, nullredMat2 = nullMat2, nullredMat3 = nullMat3,
  powerredMat2 = powerMat2, powerredMat3 = powerMat3,
  SampleSize = SampleSize, dimensions = dimensions)
}, SampleSize = SampleSize, dimensions = dimensions, reductionMethod = reductionMethod)

statistics <- filter(statistics, !(is.na(`Power Statistic`)), !(is.na(`Null Statistic`)))

save(statistics, file = paste0("E:/Ben/Box Sync/Statistics/Dissertation/Sims/", Structure, "/DimReduce/", dimensions, " ", SampleSize, "/statistics.RData"))

pushover(message = "statistics",
         title = Structure)

cvs <- summarize(group_by(statistics, SampleSize, dimension, Pops, reduction, Test),
                      CriticalValue = quantile(`Null Statistic`, 0.95))

save(cvs, file = paste0("E:/Ben/Box Sync/Statistics/Dissertation/Sims/", Structure, "/DimReduce/", dimensions, " ", SampleSize, "/cvs.RData"))

pushover(message = "cvs",
         title = Structure)


powerscorestests <- mutate(full_join(cvs, statistics),
                                Significant = (`Power Statistic` > CriticalValue))

save(powerscorestests, file = paste0("E:/Ben/Box Sync/Statistics/Dissertation/Sims/", Structure, "/DimReduce/", dimensions, " ", SampleSize, "/powerscorestests.RData"))

pushover(message = "powerscorestests",
         title = Structure)

power <- summarise(group_by(powerscorestests,
                                     SampleSize, dimension, Pops, reduction, Test),
                            Power = mean(Significant))

save(power, file = paste0("E:/Ben/Box Sync/Statistics/Dissertation/Sims/", Structure, "/DimReduce/", dimensions, " ", SampleSize, "/power.RData"))

pushover(message = "power",
         title = Structure)
