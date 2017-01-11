
power10 <- read_csv("~/Box Sync/Dissert Data/power/10.csv")
power10 <- frobenius_norm(power10)

power15 <- read_csv("~/Box Sync/Dissert Data/power/15.csv")
power15 <- frobenius_norm(power15)

power20 <- read_csv("~/Box Sync/Dissert Data/power/20.csv")
power20 <- frobenius_norm(power20)

power45 <- read_csv("~/Box Sync/Dissert Data/power/45.csv")
power45 <- frobenius_norm(power45)

power75 <- read_csv("~/Box Sync/Dissert Data/power/75.csv")
power75 <- frobenius_norm(power75)

powerboxesM <- read_csv("~/Box Sync/Dissert Data/power/boxesM.csv")
powerboxesM <- frobenius_norm(powerboxesM)

powerellip <- read_csv("~/Box Sync/Dissert Data/power/bigsimEllipticalAdd.csv")
powerellip <- frobenius_norm2(powerellip)

powertoeplitz <- read_csv("~/Box Sync/Dissert Data/power/bigsimtoeplitzmultiply.csv")
powertoeplitz <- frobenius_norm2(powertoeplitz)

power <- bind_rows(power10, power15, power20, power45, power75, powerboxesM)
power[power$type == "toeplitz",]$type <- "toeplitz add"

power <- bind_rows(power, powertoeplitz, powerellip)

reductions <- data.frame(ReductionMethod = c("DataConcatScatter",
                                             "DataConcatScatterBlock",
                                             "SDiff", "SConcat", "None"),
                         `Reduction Method` = c("Data Scatter", "Data Scatter Block",
                                                "Covariance Differences", "Covariances",
                                                "None"))

tests <- data.frame(test = c("BoxesM_test", "Chaipitak2013_test",
                             "Ishii2016_test", "Schott2007_test",
                             "Srivastava2007_test", "Srivastava2014_test",
                             "SrivastavaYanagihara2010_test"),
                    Test = c("Modified Likelihood Ratio", "Chaipitak and Chongcharoen 2013",
                             "Ishii et al. 2016", "Schott 2007",
                             "Srivastava 2007", "Srivastava et al. 2014",
                             "Srivastava and Yanagihara 2010"))
power1 <- select(full_join(full_join(power, reductions), tests), -ReductionMethod, -test)


write_csv(power1, "~/Box Sync/Dissert Data/power/powertotal.csv")
