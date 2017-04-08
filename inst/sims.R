source(system.file("simfuncEqualcov.R", package = "EqualCov"))

# sim1 <- sim(covarianceMat = toeplitz(0.5 ^ seq(0, (240 - 1))),
#             replicationscrit = 10000,
#             replicationspower = 1000,
#             samples = c(20, 40, 80, 160, 240),
#             dimensions = c(20, 40, 80, 160, 240),
#             differences = c(1.05, 1.1),
#             directory = "/Users/ben_barnard/Dropbox/toeplitz_05",
#             fileindex = "1_1a1_05")


delta0 <- matrix(0, 240, 240)

  for(i in 1:240){
    for(j in 1:240){
      delta0[i, j] <- ((-1) ^ (i + j)) * (0.4 ^ (abs(i - j) ^ 0.1))
    }
  }

omega <- diag(runif(240, 1, 5))

covmat <- omega %*% delta0 %*% omega

sim2 <- sim(covarianceMat = covmat,
            replicationscrit = 10000,
            replicationspower = 1000,
            samples = c(20, 40, 80, 160, 240),
            dimensions = c(20, 40, 80, 160, 240),
            differences = c(1.05, 1.1),
            directory = "E:/Ben/Box Sync/Statistics/unstructured/",
            fileindex = "1_1a1_05")

sim3 <- sim(covarianceMat = diag(1, 240),
            replicationscrit = 10000,
            replicationspower = 1000,
            samples = c(20, 40, 80, 160, 240),
            dimensions = c(20, 40, 80, 160, 240),
            differences = c(1.05, 1.1),
            directory = "E:/Ben/Box Sync/Statistics/identity/",
            fileindex = "1_1a1_05")
