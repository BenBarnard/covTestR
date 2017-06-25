library(MASS)
library(plyr)
library(EqualCov)
library(dplyr)
library(pushoverr)
set_pushover_user(user = "ufmfa6vc9s2fc2eh6phop9ej5ebxum")
set_pushover_app(token = "azrd3hwrwgh2gs6igbvb8yy4mftoi7")

dimensions <- c(20, 40, 80, 160)
SampleSize <- c(5, 10, 20, 40)

gridcomb <- expand.grid(Samples = SampleSize, dims = dimensions)



pushover(message = "mvndata",
         title = "Hey")

NullTests <- ldply(Map(function(dataMV){
  ldply(mapply(function(Samples, dims, df){
    listData <- list(df$Zero1[1:Samples, 1:dims], df$Zero2[1:Samples, 1:dims], df$Zero3[1:Samples, 1:dims])

    dfTest <- data.frame(SampleSize = Samples, dimension = dims,
                         Test = c("Chaipitak", "Schott", "Srivastava 2007", "Srivastava 2010",
                                  "Srivastava 2014", "Ishii", "Ahmad", "Chaipitak", "Schott",
                                  "Srivastava 2007", "Srivastava 2010", "Srivastava 2014",
                                  "Ishii", "Ahmad"),
                         Pops = c(rep("Three", 7), rep("Two", 7)),
                         Statistic = c(EqualCov:::Chaipitak2013Stat(listData),
                                       EqualCov:::Schott2007Stat(listData),
                                       EqualCov:::Srivastava2007Stat(listData),
                                       EqualCov:::SrivastavaYanagihara2010Stat(listData),
                                       EqualCov:::Srivastava2014Stat(listData),
                                       EqualCov:::Ishii2016Stat(listData),
                                       abs(EqualCov:::Ahmad2017Stat(listData)),
                                       EqualCov:::Chaipitak2013Stat(listData[-3]),
                                       EqualCov:::Schott2007Stat(listData[-3]),
                                       EqualCov:::Srivastava2007Stat(listData[-3]),
                                       EqualCov:::SrivastavaYanagihara2010Stat(listData[-3]),
                                       EqualCov:::Srivastava2014Stat(listData[-3]),
                                       EqualCov:::Ishii2016Stat(listData[-3]),
                                       abs(EqualCov:::Ahmad2017Stat(listData[-3]))))
    dfTest
  }, Samples = gridcomb$Samples, dims = gridcomb$dims, MoreArgs = list(df = dataMV), SIMPLIFY = FALSE))
}, mvndata))

save(NullTests, file = "E:/Ben/Box Sync/Statistics/Autoregressive/NullTests.RData")

pushover(message = "NullTests",
         title = "Hey")

cvs <- summarize(group_by(NullTests, SampleSize, dimension, Test, Pops),
                 CriticalValue = quantile(Statistic, 0.95))

save(cvs, file = "E:/Ben/Box Sync/Statistics/Autoregressive/cvs.RData")

pushover(message = "cvs",
         title = "Hey")

PowerTests <- ldply(Map(function(dataMV){
  ldply(mapply(function(Samples, dims, df){
    listData <- list(df$Zero1[1:Samples, 1:dims], df$One[1:Samples, 1:dims], df$Two[1:Samples, 1:dims])

    dfTest <- data.frame(SampleSize = Samples, dimension = dims,
                         Test = c("Chaipitak", "Schott", "Srivastava 2007", "Srivastava 2010",
                                  "Srivastava 2014", "Ishii", "Ahmad", "Chaipitak", "Schott",
                                  "Srivastava 2007", "Srivastava 2010", "Srivastava 2014",
                                  "Ishii", "Ahmad"),
                         Pops = c(rep("Three", 7), rep("Two", 7)),
                         Statistic = c(EqualCov:::Chaipitak2013Stat(listData),
                                       EqualCov:::Schott2007Stat(listData),
                                       EqualCov:::Srivastava2007Stat(listData),
                                       EqualCov:::SrivastavaYanagihara2010Stat(listData),
                                       EqualCov:::Srivastava2014Stat(listData),
                                       EqualCov:::Ishii2016Stat(listData),
                                       abs(EqualCov:::Ahmad2017Stat(listData)),
                                       EqualCov:::Chaipitak2013Stat(listData[-3]),
                                       EqualCov:::Schott2007Stat(listData[-3]),
                                       EqualCov:::Srivastava2007Stat(listData[-3]),
                                       EqualCov:::SrivastavaYanagihara2010Stat(listData[-3]),
                                       EqualCov:::Srivastava2014Stat(listData[-3]),
                                       EqualCov:::Ishii2016Stat(listData[-3]),
                                       abs(EqualCov:::Ahmad2017Stat(listData[-3]))))
    dfTest
  }, Samples = gridcomb$Samples, dims = gridcomb$dims, MoreArgs = list(df = dataMV), SIMPLIFY = FALSE))
}, mvndata))

save(PowerTests, file = "E:/Ben/Box Sync/Statistics/Autoregressive/PowerTests.RData")

pushover(message = "PowerTests",
         title = "Hey")

powerscorestests <- mutate(full_join(cvs, PowerTests),
                           Significant = (Statistic > CriticalValue))

save(powerscorestests, file = "E:/Ben/Box Sync/Statistics/Autoregressive/powerscorestests.RData")

pushover(message = "powerscorestests",
         title = "Hey")

power <- summarise(group_by(powerscorestests,
                            SampleSize, dimension, Pops, Test),
                   Power = mean(Significant))

save(power, file = "E:/Ben/Box Sync/Statistics/Autoregressive/power.RData")

pushover(message = "power",
         title = "Hey")
