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

multi <- function(dimensions, samples, critdat){

  critdat <- critdat[names(critdat) == samples][[1]]

  data_frame(Test = c(rep("Chaipitak", length(critdat)),
                      rep("Schott", length(critdat)),
                      rep("Srivastava 2007", length(critdat)),
                      rep("Srivastava 2010", length(critdat)),
                      rep("Srivastava 2014", length(critdat)),
                      rep("Ishii", length(critdat)),
                      rep("Chaipitak", length(critdat)),
                      rep("Schott", length(critdat)),
                      rep("Srivastava 2007", length(critdat)),
                      rep("Srivastava 2010", length(critdat)),
                      rep("Srivastava 2014", length(critdat)),
                      rep("Ishii", length(critdat))),
             Value = c(ldply(map(.x = critdat, test3, test = Chaipitak2013_test, dimensions = dimensions))$V1,
                       ldply(map(.x = critdat, test3, test = Schott2007_test, dimensions = dimensions))$V1,
                       ldply(map(.x = critdat, test3, test = Srivastava2007_test, dimensions = dimensions))$V1,
                       ldply(map(.x = critdat, test3, test = SrivastavaYanagihara2010_test, dimensions = dimensions))$V1,
                       ldply(map(.x = critdat, test3, test = Srivastava2014_test, dimensions = dimensions))$V1,
                       ldply(map(.x = critdat, test3, test = Ishii2016_test, dimensions = dimensions))$V1,
                       ldply(map(.x = critdat, test2, test = Chaipitak2013_test, dimensions = dimensions))$V1,
                       ldply(map(.x = critdat, test2, test = Schott2007_test, dimensions = dimensions))$V1,
                       ldply(map(.x = critdat, test2, test = Srivastava2007_test, dimensions = dimensions))$V1,
                       ldply(map(.x = critdat, test2, test = SrivastavaYanagihara2010_test, dimensions = dimensions))$V1,
                       ldply(map(.x = critdat, test2, test = Srivastava2014_test, dimensions = dimensions))$V1,
                       ldply(map(.x = critdat, test2, test = Ishii2016_test, dimensions = dimensions))$V1),
             Populations = c(rep(3, 6 * length(critdat)), rep(2, 6 * length(critdat))))

}

powerdata <- function(samples, difference, covarianceMat, replicationspower){
  replicate(
    n = replicationspower,
    expr = c(
      wishart::rWishart(n = 1,
                        df = samples - 1,
                        Sigma = covarianceMat,
                        covariance = TRUE, simplify = FALSE),
      wishart::rWishart(n = 1,
                        df = samples - 1,
                        Sigma = covarianceMat * difference,
                        covariance = TRUE, simplify = FALSE),
      wishart::rWishart(n = 1,
                        df = samples - 1,
                        Sigma = covarianceMat,
                        covariance = TRUE, simplify = FALSE)),
    simplify = FALSE)
}

multipower <- function(dimensions, samples, differences, powerdat){


  labels <- attributes(powerdat)$split_labels

  grp <- which(labels$samples == samples & labels$difference == 1.1)

  powerdat <- powerdat[names(powerdat) == grp][[1]]

  data_frame(Test = c(rep("Chaipitak", length(powerdat)),
                      rep("Schott", length(powerdat)),
                      rep("Srivastava 2007", length(powerdat)),
                      rep("Srivastava 2010", length(powerdat)),
                      rep("Srivastava 2014", length(powerdat)),
                      rep("Ishii", length(powerdat)),
                      rep("Chaipitak", length(powerdat)),
                      rep("Schott", length(powerdat)),
                      rep("Srivastava 2007", length(powerdat)),
                      rep("Srivastava 2010", length(powerdat)),
                      rep("Srivastava 2014", length(powerdat)),
                      rep("Ishii", length(powerdat))),
             Value = c(ldply(map(.x = powerdat, test3, test = Chaipitak2013_test, dimensions = dimensions))$V1,
                       ldply(map(.x = powerdat, test3, test = Schott2007_test, dimensions = dimensions))$V1,
                       ldply(map(.x = powerdat, test3, test = Srivastava2007_test, dimensions = dimensions))$V1,
                       ldply(map(.x = powerdat, test3, test = SrivastavaYanagihara2010_test, dimensions = dimensions))$V1,
                       ldply(map(.x = powerdat, test3, test = Srivastava2014_test, dimensions = dimensions))$V1,
                       ldply(map(.x = powerdat, test3, test = Ishii2016_test, dimensions = dimensions))$V1,
                       ldply(map(.x = powerdat, test2, test = Chaipitak2013_test, dimensions = dimensions))$V1,
                       ldply(map(.x = powerdat, test2, test = Schott2007_test, dimensions = dimensions))$V1,
                       ldply(map(.x = powerdat, test2, test = Srivastava2007_test, dimensions = dimensions))$V1,
                       ldply(map(.x = powerdat, test2, test = SrivastavaYanagihara2010_test, dimensions = dimensions))$V1,
                       ldply(map(.x = powerdat, test2, test = Srivastava2014_test, dimensions = dimensions))$V1,
                       ldply(map(.x = powerdat, test2, test = Ishii2016_test, dimensions = dimensions))$V1),
             Populations = c(rep(3, 6 * length(powerdat)), rep(2, 6 * length(powerdat))))

}


powervalue = function(powertests, critValues){
  critValue <- filter(critValues, dimensions == unique(powertests$dimensions), samples == unique(powertests$samples),
                      Test == unique(powertests$Test), Populations == unique(powertests$Populations))$`Critical Value`

  powertests <- mutate(powertests, Cond = ifelse(Value >= critValue, 1, 0))

  summarise(group_by(powertests, samples, dimensions, difference, Test, Populations), Power = mean(Cond))
}

