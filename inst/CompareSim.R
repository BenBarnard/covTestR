library(MASS)
library(plyr)
library(EqualCov)
library(dplyr)
library(tibble)
library(pushoverr)
set_pushover_user(user = "ufmfa6vc9s2fc2eh6phop9ej5ebxum")
set_pushover_app(token = "azrd3hwrwgh2gs6igbvb8yy4mftoi7")

dimensions <- c(20, 40, 80, 160)
SampleSize <- c(5, 10, 20, 40)
covStruct <- "Unstructured"
replications <- 1000

gridcomb <- filter(expand.grid(Samples = SampleSize,
                               dims = dimensions),
                   dims > Samples)


load(paste0("E:/Ben/Box Sync/Statistics/Dissertation/mvnData/",
            covStruct, "/mvndata.RData"))

mvndata <- mvndata[1:replications]

#### Statistics ####

statistics <- ldply(lapply(mvndata, function(datadf, grid){
  ldply(mapply(function(Samp, Dim, df){

    nullList <- list(datadf$Zero1[1:Samp, 1:Dim],
                     datadf$Zero2[1:Samp, 1:Dim],
                     datadf$Zero3[1:Samp, 1:Dim])

    powerList <- list(datadf$Zero1[1:Samp, 1:Dim],
                      datadf$One[1:Samp, 1:Dim],
                      datadf$Two[1:Samp, 1:Dim])

    dfTest <- data_frame(SampleSize = Samp, dimension = Dim,
                         Test = rep(c("Chaipitak", "Chaipitak pool", "Schott", "Schott pool",
                                      "Srivastava 2007", "Srivastava 2010", "Srivastava 2014",
                                      "Srivastava 2014 pool", "Ishii", "Ahmad"), 2),
                         Pops = c(rep("Three", 10), rep("Two", 10)),
                         `Null Statistic` = c(EqualCov:::Chaipitak2013Stat(nullList),
                                              EqualCov:::Chaipitak2013poolStat(nullList),
                                              EqualCov:::Schott2007Stat(nullList),
                                              EqualCov:::Schott2007pooledStat(nullList),
                                              EqualCov:::Srivastava2007Stat(nullList),
                                              EqualCov:::SrivastavaYanagihara2010Stat(nullList),
                                              EqualCov:::Srivastava2014Stat(nullList),
                                              EqualCov:::Srivastava2014poolStat(nullList),
                                              EqualCov:::Ishii2016Stat(nullList),
                                              abs(EqualCov:::Ahmad2017Stat(nullList)),
                                              EqualCov:::Chaipitak2013Stat(nullList[-3]),
                                              EqualCov:::Chaipitak2013poolStat(nullList[-3]),
                                              EqualCov:::Schott2007Stat(nullList[-3]),
                                              EqualCov:::Schott2007pooledStat(nullList[-3]),
                                              EqualCov:::Srivastava2007Stat(nullList[-3]),
                                              EqualCov:::SrivastavaYanagihara2010Stat(nullList[-3]),
                                              EqualCov:::Srivastava2014Stat(nullList[-3]),
                                              EqualCov:::Srivastava2014poolStat(nullList[-3]),
                                              EqualCov:::Ishii2016Stat(nullList[-3]),
                                              abs(EqualCov:::Ahmad2017Stat(nullList[-3]))),
                         `Power Statistic` = c(EqualCov:::Chaipitak2013Stat(powerList),
                                               EqualCov:::Chaipitak2013poolStat(powerList),
                                               EqualCov:::Schott2007Stat(powerList),
                                               EqualCov:::Schott2007pooledStat(powerList),
                                               EqualCov:::Srivastava2007Stat(powerList),
                                               EqualCov:::SrivastavaYanagihara2010Stat(powerList),
                                               EqualCov:::Srivastava2014Stat(powerList),
                                               EqualCov:::Srivastava2014poolStat(powerList),
                                               EqualCov:::Ishii2016Stat(powerList),
                                               abs(EqualCov:::Ahmad2017Stat(powerList)),
                                               EqualCov:::Chaipitak2013Stat(powerList[-3]),
                                               EqualCov:::Chaipitak2013poolStat(powerList[-3]),
                                               EqualCov:::Schott2007Stat(powerList[-3]),
                                               EqualCov:::Schott2007pooledStat(powerList[-3]),
                                               EqualCov:::Srivastava2007Stat(powerList[-3]),
                                               EqualCov:::SrivastavaYanagihara2010Stat(powerList[-3]),
                                               EqualCov:::Srivastava2014Stat(powerList[-3]),
                                               EqualCov:::Srivastava2014poolStat(powerList[-3]),
                                               EqualCov:::Ishii2016Stat(powerList[-3]),
                                               abs(EqualCov:::Ahmad2017Stat(powerList[-3]))))

    dfTest
  }, Samp = grid$Samples, Dim = grid$dims,
  MoreArgs = list(df = datadf), SIMPLIFY = FALSE))
}, grid = gridcomb))


save(statistics, file = paste0("E:/Ben/Box Sync/Statistics/Dissertation/",
                               covStruct, "/statistics.RData"))

pushover(message = "Statistics",
         title = covStruct)


#### Critical Values ####

cvs <- summarize(group_by(statistics, SampleSize, dimension, Test, Pops),
                 CriticalValue = quantile(`Null Statistic`, 0.95))

save(cvs, file = paste0("E:/Ben/Box Sync/Statistics/Dissertation/",
                        covStruct, "/cvs.RData"))

pushover(message = "cvs",
         title = covStruct)



#### Power Score Tests ####

powerscorestests <- mutate(full_join(cvs, statistics),
                           Significant = (`Power Statistic` > CriticalValue))

save(powerscorestests, file = paste0("E:/Ben/Box Sync/Statistics/Dissertation/",
                                     covStruct, "/powerscoretests.RData"))

pushover(message = "powerscorestests",
         title = covStruct)


#### Power ####

power <- summarise(group_by(powerscorestests,
                            SampleSize, dimension, Pops, Test),
                   Power = mean(Significant))

save(power, file = paste0("E:/Ben/Box Sync/Statistics/Dissertation/",
                          covStruct, "/power.RData"))

pushover(message = "power",
         title = covStruct)
