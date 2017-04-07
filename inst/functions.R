library(plyr)
library(dplyr)
library(EqualCov)
library(MASS)
library(purrr)

test3 <- function(data, test, dimensions){
  data <- lapply(data, function(x){
    x[1:dimensions, 1:dimensions]
  })
  do.call(test, list(data))$statistic[[1]]
}

test2 <- function(data, test, dimensions){
  data <- lapply(data, function(x){
    x[1:dimensions, 1:dimensions]
  })
  do.call(test, list(data[1:2]))$statistic[[1]]
}

multi <- function(x, dimensions, covarianceMat){
    dat <- wishart::rWishart(
      n = 3,
      df = x - 1,
      Sigma = covarianceMat,
      covariance = TRUE, simplify = FALSE)

    dimensions <- dimensions[dimensions >= x]

ldply(dimensions, function(dimensions, dat, samples){
  data_frame(Test =c("Chaipitak", "Schott", "Srivastava 2007", "Srivastava 2010",
                     "Srivastava 2014", "Ishii", "Chaipitak", "Schott",
                     "Srivastava 2007", "Srivastava 2010", "Srivastava 2014", "Ishii"),
             Value = c(
               test3(dat, test = Chaipitak2013_test, dimensions = dimensions),
               test3(dat, test = Schott2007_test, dimensions = dimensions),
               test3(dat, test = Srivastava2007_test, dimensions = dimensions),
               test3(dat, test = SrivastavaYanagihara2010_test, dimensions = dimensions),
               test3(dat, test = Srivastava2014_test, dimensions = dimensions),
               test3(dat, test = Ishii2016_test, dimensions = dimensions),
               test2(dat, test = Chaipitak2013_test, dimensions = dimensions),
               test2(dat, test = Schott2007_test, dimensions = dimensions),
               test2(dat, test = Srivastava2007_test, dimensions = dimensions),
               test2(dat, test = SrivastavaYanagihara2010_test, dimensions = dimensions),
               test2(dat, test = Srivastava2014_test, dimensions = dimensions),
               test2(dat, test = Ishii2016_test, dimensions = dimensions)),
             Populations = c(rep(3, 6), rep(2, 6)),
             Samples = rep(samples, 12),
             Dimensions = dimensions)}, dat = dat, samples = x,
  .progress = "text")
}

powerdata <- function(dimensions, Samples, Differences, covarianceMat, replicationspower){
  ldply(replicate(
    n = replicationspower,
    expr = multipower(dimensions, Samples, Differences, covarianceMat),
    simplify = FALSE))
}

multipower <- function(dimensions, samples, differences, covarianceMat){

  dat <- c(wishart::rWishart(n = 1,
                      df = samples - 1,
                      Sigma = covarianceMat,
                      covariance = TRUE, simplify = FALSE),
    wishart::rWishart(n = 1,
                      df = samples - 1,
                      Sigma = covarianceMat * differences,
                      covariance = TRUE, simplify = FALSE),
    wishart::rWishart(n = 1,
                      df = samples - 1,
                      Sigma = covarianceMat,
                      covariance = TRUE, simplify = FALSE))

  dimensions <- dimensions[dimensions >= samples]

  ldply(dimensions, function(dimensions, dat, samples, differences){
    data_frame(Test =c("Chaipitak", "Schott", "Srivastava 2007", "Srivastava 2010",
                       "Srivastava 2014", "Ishii", "Chaipitak", "Schott",
                       "Srivastava 2007", "Srivastava 2010", "Srivastava 2014", "Ishii"),
               Value = c(
                 test3(dat, test = Chaipitak2013_test, dimensions = dimensions),
                 test3(dat, test = Schott2007_test, dimensions = dimensions),
                 test3(dat, test = Srivastava2007_test, dimensions = dimensions),
                 test3(dat, test = SrivastavaYanagihara2010_test, dimensions = dimensions),
                 test3(dat, test = Srivastava2014_test, dimensions = dimensions),
                 test3(dat, test = Ishii2016_test, dimensions = dimensions),
                 test2(dat, test = Chaipitak2013_test, dimensions = dimensions),
                 test2(dat, test = Schott2007_test, dimensions = dimensions),
                 test2(dat, test = Srivastava2007_test, dimensions = dimensions),
                 test2(dat, test = SrivastavaYanagihara2010_test, dimensions = dimensions),
                 test2(dat, test = Srivastava2014_test, dimensions = dimensions),
                 test2(dat, test = Ishii2016_test, dimensions = dimensions)),
               Populations = c(rep(3, 6), rep(2, 6)),
               Dimensions = dimensions)}, dat = dat, samples = samples, differences, .progress = progress_text(char = "q"))
}


powervalue = function(powertests, critValues){
  critValue <- filter(critValues, Dimensions == unique(powertests$Dimensions), Samples == unique(powertests$Samples),
                      Test == unique(powertests$Test), Populations == unique(powertests$Populations))$`Critical Value`

  powertests <- mutate(powertests, Cond = ifelse(Value >= critValue, 1, 0))

  summarise(group_by(powertests, Samples, Dimensions, Differences, Test, Populations), Power = mean(Cond))
}

