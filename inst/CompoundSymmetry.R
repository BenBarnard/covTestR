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

Sigma <- 0.4 * diag(1, 200) + 0.6 * rep(1, 200) %*% t(rep(1, 200))

mvndata <- replicate(replications, list("CS" = mvrnorm(n = SampleSize, mu = rep(0, dimensions),
                                                 Sigma = Sigma),
                               "Ident" = mvrnorm(n = SampleSize, mu = rep(0, dimensions),
                                                 Sigma = diag(1, dimensions))),
            simplify = FALSE)

pushover(message = "mvndata",
         title = "Hey")

NullTests <- ldply(mvndata, function(ls, SampleSize, dimensions){

  redMat <- songStructure(ls$Ident)$u
  ldply(c(4:SampleSize, dimensions), function(reduction, ls, redMat){

    lt <- ls
    SampleSize <- nrow(lt$CS)
    originaldimension <- ncol(lt$CS)
    projection <- t(redMat[,1:reduction])

    if(reduction == originaldimension){
      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                       reduction = reduction,
                       Test = c("Ahmad 2015", "Chen 2010", "Fisher 2012", "Ledoit 2002",
                                "Nagao 1973", "Srivastava 2005", "Srivastava 2011"),
                       Statistic = c(abs(Ahmad2015(lt$Ident)$statistic),
                                     abs(Chen2010(lt$Ident)$statistic),
                                     abs(Fisher2012(lt$Ident)$statistic),
                                     abs(LedoitWolf2002(lt$Ident)$statistic),
                                     abs(Nagao1973(lt$Ident)$statistic),
                                     abs(Srivastava2005(lt$Ident)$statistic),
                                     abs(Srivastava2011(lt$Ident)$statistic)))

    }

    lt <- lapply(lt, function(x){
      t(projection %*% t(x))
    })

    proj <- cov(t(projection))

    if(reduction <= SampleSize){
      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                       reduction = reduction,
                       Test = c("Ahmad 2015", "Chen 2010", "Fisher 2012", "Ledoit 2002",
                                "Nagao 1973", "Srivastava 2005", "Srivastava 2011"),
                       Statistic = c(abs(Ahmad2015(lt$Ident, Sigma = proj)$statistic),
                                     abs(Chen2010(lt$Ident, Sigma = proj)$statistic),
                                     abs(Fisher2012(lt$Ident, Sigma = proj)$statistic),
                                     abs(LedoitWolf2002(lt$Ident, Sigma = proj)$statistic),
                                     abs(Nagao1973(lt$Ident, Sigma = proj)$statistic),
                                     abs(Srivastava2005(lt$Ident, Sigma = proj)$statistic),
                                     abs(Srivastava2011(lt$Ident, Sigma = proj)$statistic)))

    }

    df

  }, ls = ls, redMat = redMat)
}, SampleSize = SampleSize, dimensions = dimensions)

save(NullTests, file = "E:/Ben/Box Sync/Statistics/OneSampleSimDimReduce/NullTests.RData")

pushover(message = "NullTests",
         title = "Hey")

cvs <- summarize(group_by(NullTests, SampleSize, originaldimension, reduction, Test),
                      CriticalValue = quantile(Statistic, 0.95))

save(cvs, file = "E:/Ben/Box Sync/Statistics/OneSampleSimDimReduce/cvs.RData")

pushover(message = "cvs",
         title = "Hey")

powerValueTests <- ldply(mvndata, function(ls, SampleSize, dimensions){

  redMat <- songStructure(ls$CS)$u
  ldply(c(4:SampleSize, dimensions), function(reduction, ls, redMat){

    lt <- ls
    SampleSize <- nrow(lt$CS)
    originaldimension <- ncol(lt$CS)
    projection <- t(redMat[,1:reduction])

    if(reduction == originaldimension){
      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                       reduction = reduction,
                       Test = c("Ahmad 2015", "Chen 2010", "Fisher 2012", "Ledoit 2002",
                                "Nagao 1973", "Srivastava 2005", "Srivastava 2011"),
                       Statistic = c(abs(Ahmad2015(lt$CS)$statistic),
                                     abs(Chen2010(lt$CS)$statistic),
                                     abs(Fisher2012(lt$CS)$statistic),
                                     abs(LedoitWolf2002(lt$CS)$statistic),
                                     abs(Nagao1973(lt$CS)$statistic),
                                     abs(Srivastava2005(lt$CS)$statistic),
                                     abs(Srivastava2011(lt$CS)$statistic)))

    }

    lt <- lapply(lt, function(x){
      t(projection %*% t(x))
    })

    proj <- cov(t(projection))

    if(reduction <= SampleSize){
      df <- data.frame(SampleSize = SampleSize, originaldimension = originaldimension,
                       reduction = reduction,
                       Test = c("Ahmad 2015", "Chen 2010", "Fisher 2012", "Ledoit 2002",
                                "Nagao 1973", "Srivastava 2005", "Srivastava 2011"),
                       Statistic = c(abs(Ahmad2015(lt$CS, Sigma = proj)$statistic),
                                     abs(Chen2010(lt$CS, Sigma = proj)$statistic),
                                     abs(Fisher2012(lt$CS, Sigma = proj)$statistic),
                                     abs(LedoitWolf2002(lt$CS, Sigma = proj)$statistic),
                                     abs(Nagao1973(lt$CS, Sigma = proj)$statistic),
                                     abs(Srivastava2005(lt$CS, Sigma = proj)$statistic),
                                     abs(Srivastava2011(lt$CS, Sigma = proj)$statistic)))

    }

    df

  }, ls = ls, redMat = redMat)
}, SampleSize = SampleSize, dimensions = dimensions)

save(powerValueTests, file = "E:/Ben/Box Sync/Statistics/OneSampleSimDimReduce/powerValueTests.RData")

pushover(message = "powerValueTests",
         title = "Hey")


powerscorestests <- mutate(full_join(cvs, powerValueTests),
                                Significant = (Statistic > CriticalValue))

save(powerscorestests, file = "E:/Ben/Box Sync/Statistics/OneSampleSimDimReduce/powerscorestests.RData")

pushover(message = "powerscorestests",
         title = "Hey")

powertest <- summarise(group_by(powerscorestests,
                                     SampleSize, originaldimension, reduction, Test),
                            Power = mean(Significant))

save(powerthreetest, file = "E:/Ben/Box Sync/Statistics/OneSampleSimDimReduce/powertest.RData")

pushover(message = "powertest",
         title = "Hey")

