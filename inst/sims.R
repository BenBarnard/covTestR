source(system.file("simfuncEqualcov.R", package = "EqualCov"))

sim1 <- sim(covarianceMat = diag(1, 100),
            replicationscrit = 10,
            repilcationspower = 5,
            samples = c(20, 30),
            dimensions = c(80, 100),
            differences = 1.05,
            directory = "/Users/ben_barnard/Dropbox/")
