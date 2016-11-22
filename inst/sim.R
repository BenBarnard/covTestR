crit_data(type = "identity",
          dimensions = c(100),
          samples = c(15, 45, 75),
          difference = seq(.1, .9, .1),
          populations = 2,
          replications = 1000,
          reductionMethods = c("DataConcatScatter", "DataConcatScatterBlock",
                               "SConcat", "SDiff"),
          directory = "~/Documents/R/Dissertation/data/100")

crit_data3(type = "identity",
           dimensions = c(100),
           samples = c(15, 45, 75),
           difference = seq(.1, .9, .1),
           populations = 3,
           replications = 1000,
           reductionMethods = c("DataConcatScatter", "DataConcatScatterBlock",
                                "SConcat", "SDiff"),
           directory = "~/Documents/R/Dissertation/data/100")

crit_data(type = "toeplitz",
          dimensions = c(100),
          samples = c(15, 45, 75),
          difference = seq(.1, .9, .1),
          populations = 2,
          replications = 1000,
          reductionMethods = c("DataConcatScatter", "DataConcatScatterBlock",
                               "SConcat", "SDiff"),
          directory = "~/Documents/R/Dissertation/data/100")

crit_data3(type = "toeplitz",
           dimensions = c(100),
           samples = c(15, 45, 75),
           difference = seq(.1, .9, .1),
           populations = 3,
           replications = 1000,
           reductionMethods = c("DataConcatScatter", "DataConcatScatterBlock",
                                "SConcat", "SDiff"),
           directory = "~/Documents/R/Dissertation/data/100")

tests(directory = "~/Documents/R/Dissertation/data/100",
      save = "~/Dropbox",
      test_funcs = c("Schott2007_test", "Srivastava2007_test",
                     "SrivastavaYanagihara2010_test",
                     "Srivastava2014_test", "Ishii2016_test",
                     "Chaipitak2013_test"),
      dimensions = seq(20, 5))
