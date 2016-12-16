crit_data(type = "toeplitz ",
          dimensions = c(100, 50),
          samples = c(10, 15, 20),
          difference = seq(0, 2, .1),
          populations = 2,
          replications = 1000,
          reductionMethods = c("DataConcatScatter",
                               "SConcat", "SDiff"),
          directory = "~/Documents/R/Dissertation/data/100_5")

crit_data3(type = "toeplitz ",
           dimensions = c(100, 50),
           samples = c(10, 15, 20),
           difference = seq(0, 2, .1),
           populations = 3,
           replications = 1000,
           reductionMethods = c("DataConcatScatter",
                                "SConcat", "SDiff"),
           directory = "~/Documents/R/Dissertation/data/100_5")

tests(directory = "~/Documents/R/Dissertation/data/100_5",
      save = "~/Dropbox",
      test_funcs = c("Schott2007_test", "Srivastava2007_test",
                     "SrivastavaYanagihara2010_test",
                     "Srivastava2014_test", "Ishii2016_test",
                     "Chaipitak2013_test", "BoxesM_test"),
      dimensions = seq(20, 5))

