sim <- function(covarianceMat, replicationscrit, replicationspower, samples, dimensions, differences, directory){

  source(system.file("functions.R", package = "EqualCov"))
  browser()
  samplescrit <- setNames(samples, samples)

  critical_data <- llply(samplescrit,
                         function(x, covarianceMat, replicationscrit){
                           replicate(
                             n = replicationscrit,
                             expr = wishart::rWishart(
                               n = 3,
                               df = x - 1,
                               Sigma = covarianceMat,
                               covariance = TRUE, simplify = FALSE), simplify = FALSE
                           )
                         }, covarianceMat = covarianceMat,
                         replicationscrit = replicationscrit,
                         .progress = "text")

  save(critical_data, file = paste0(directory, "critical_data.RData"))

  critgrid <- filter(expand.grid(samples = samples, dimensions = dimensions), samples <= dimensions)

  critical_tests <- mdply(.data = critgrid,
                          .fun = multi, .progress = progress_text(char = "p"), critdat = critical_data)

  save(critical_tests, file = paste0(directory, "critical_tests.RData"))

  critical_values <- summarize(group_by(critical_tests, dimensions, samples, Test, Populations), `Critical Value` = quantile(Value, 0.95))

  save(critical_values, file = "critical_values.RData")

  powerdatgrid <- expand.grid(samples = samples, difference = differences)

  power_data <- mlply(powerdatgrid,
                      powerdata, covarianceMat = covarianceMat, replicationspower = replicationspower, .progress = progress_text(char = "-"))

  save(power_data, file = paste0(directory, "power_data.RData"))

  powertestgrid <- filter(expand.grid(samples = samples, dimensions = dimensions, difference = differences), samples <= dimensions)

  power_tests <- mdply(powertestgrid,
                       .fun = multipower, .progress = progress_text(char = "q"), powerdat = power_data)

  save(power_tests, file = paste0(directory, "power_tests.RData"))

  power_values <- ddply(.data = power_tests, .variables = .(samples, dimensions, difference, Test, Populations),
                        .fun = powervalue, .progress = progress_text(char = "u"), critValues = critical_values)

  save(power_values, file = paste0(directory, "power_values.RData"))

}
