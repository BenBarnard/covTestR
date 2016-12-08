crit_data(type = "toeplitz",
          dimensions = c(100),
          samples = c(10, 20),
          difference = 0,
          populations = 2,
          replications = 1000,
          reductionMethods = c("DataConcatScatter", "DataConcatScatterBlock",
                               "SConcat", "SDiff"),
          directory = "~/Documents/R/Dissertation/data/100_3")

crit_data3(type = "toeplitz",
           dimensions = c(100),
           samples = c(10, 20),
           difference = 0,
           populations = 3,
           replications = 1000,
           reductionMethods = c("DataConcatScatter", "DataConcatScatterBlock",
                                "SConcat", "SDiff"),
           directory = "~/Documents/R/Dissertation/data/100_3")

crit_data(type = "identity",
          dimensions = c(100),
          samples = c(10, 20),
          difference = 0,
          populations = 2,
          replications = 1000,
          reductionMethods = c("DataConcatScatter", "DataConcatScatterBlock",
                               "SConcat", "SDiff"),
          directory = "~/Documents/R/Dissertation/data/100_3")

crit_data3(type = "identity",
           dimensions = c(100),
           samples = c(10, 20),
           difference = 0,
           populations = 3,
           replications = 1000,
           reductionMethods = c("DataConcatScatter", "DataConcatScatterBlock",
                                "SConcat", "SDiff"),
           directory = "~/Documents/R/Dissertation/data/100_3")

tests(directory = "~/Documents/R/Dissertation/data/100_3",
      save = "~/Dropbox",
      test_funcs = c("Schott2007_test", "Srivastava2007_test",
                     "SrivastavaYanagihara2010_test",
                     "Srivastava2014_test", "Ishii2016_test",
                     "Chaipitak2013_test", "BoxesM_test"),
      dimensions = seq(20, 5))
