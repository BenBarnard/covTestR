sim <- function(covarianceMat, replicationscrit, replicationspower, samples, dimensions, differences, directory){

  source(system.file("functions.R", package = "EqualCov"))

files <- list.files(directory)

if(!("critical_tests.RData" %in% files)){
  critical_tests <- ldply(samples,
                         function(x, covarianceMat, replicationscrit, dimensions){
                           ldply(replicate(
                             n = replicationscrit,
                             expr = multi(x, dimensions, covarianceMat),
                             simplify = FALSE))
                         }, covarianceMat = covarianceMat,
                         replicationscrit = replicationscrit,
                         dimensions = dimensions)
  save(critical_tests, file = paste0(directory, "critical_tests.RData"))
}else{
  load(paste0(directory, "critical_tests.RData"))
}


if(!("critical_values.RData" %in% files)){
  critical_values <- summarize(group_by(critical_tests, Dimensions, Samples, Test, Populations), `Critical Value` = quantile(Value, 0.95))

  save(critical_values, file = paste0(directory, "critical_values.RData"))
}else{
  load(paste0(directory, "critical_values.RData"))
}

powertestgrid <- expand.grid(Samples = samples, Differences = differences)


if(!("power_tests.RData" %in% files)){

  power_tests <- mdply(powertestgrid,
                       .fun = powerdata,
                       covarianceMat = covarianceMat, replicationspower = replicationspower, dimensions = dimensions)

  save(power_tests, file = paste0(directory, "power_tests", fileindex, ".RData"))
}else{
  load(paste0(directory, "power_tests.RData"))
browser()
}

if(!("power_tests.RData" %in% files)){
  power_values <- ddply(.data = power_tests, .variables = .(Samples, Dimensions, Differences, Test, Populations),
                        .fun = powervalue, .progress = progress_text(char = "u"), critValues = critical_values)

  save(power_values, file = paste0(directory, "power_values", fileindex, ".RData"))
}else{
  load(paste0(directory, "power_tests.RData"))
  browser()
}

  power_values
}
