sim <- function(covarianceMat, replicationscrit, repilcationspower, samples, dimensions, differences, directory){

setwd(directory)
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

save(critical_data, file = "critical_data.RData")

critgrid <- filter(expand.grid(samples = samples, dimensions = dimensions), samples <= dimensions)

critical_tests <- mdply(.data = critgrid,
                        .fun = multi, .progress = progress_text(char = "p"), critdat = critical_data)

save(critical_tests, file = "critical_tests.RData")

critical_values <- summarize(group_by(critical_tests, dimensions, samples, Test, Populations), `Critical Value` = quantile(Value, 0.95))

save(critical_values, file = "critical_values.RData")

power_data <- mlply(data_frame(samples = c(20, 40, 80, 160, 240, 20, 40, 80, 160, 240),
                               difference = c(rep(1.05, 5), rep(1.1, 5))),
                    powerdata, covarianceMat = covarianceMat, replicationspower = replicationspower, .progress = progress_text(char = "-"))

save(power_data, file = "power_data.RData")

power_tests <- mdply(data_frame(samples = c(20, 20, 40, 20, 40, 80, 20, 40, 80, 160, 20, 40, 80, 160, 240,
                                            20, 20, 40, 20, 40, 80, 20, 40, 80, 160, 20, 40, 80, 160, 240),
                                dimensions = c(20, 40, 40, 80, 80, 80, 160, 160, 160, 160, 240, 240, 240, 240, 240,
                                               20, 40, 40, 80, 80, 80, 160, 160, 160, 160, 240, 240, 240, 240, 240),
                                difference = c(rep(1.05, 15), rep(1.1, 15))),
                     .fun = multipower, .progress = progress_text(char = "q"), powerdat = power_data)

save(power_tests, file = "power_tests.RData")

power_values <- ddply(.data = power_tests, .variables = .(samples, dimensions, difference, Test, Populations),
                      .fun = powervalue, .progress = progress_text(char = "u"), critValues = critical_values)

save(power_values, file = "power_values.RData")

}
