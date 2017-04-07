source(system.file("simfuncEqualcov.R", package = "EqualCov"))

sim1 <- sim(covarianceMat = toeplitz(0.5 ^ seq(0, (240 - 1))),
            replicationscrit = 10000,
            replicationspower = 1000,
            samples = c(20, 40, 80, 160, 240),
            dimensions = c(20, 40, 80, 160, 240),
            differences = c(1.05, 1.1),
            directory = "/Users/ben_barnard/Dropbox/toeplitz_05",
            fileindex = "1_1a1_05")
