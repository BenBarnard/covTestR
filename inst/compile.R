library(plyr)
library(stringr)
library(tidyverse)

files <- list.files("~/Box Sync/Dissert Data/100")

data <- ldply(files, function(x){
  filegroup <- as.data.frame(str_split(x, " ", simplify = TRUE))
  names(filegroup) <- c("dimension", "samples", "difference", "type", "extension")
  data <- read_csv(paste0("~/Box Sync/Dissert Data/100/", x))
  cbind(data, data.frame(Samples = rep(filegroup$samples, nrow(data))))
})

critical_data <- filter(data, difference == 0)

critical_values <- summarise(group_by(critical_data, ReductionMethod, type,
                                      ReducedDimension, populations, Samples, test),
                             Crit = quantile(value, probs = .95))

power_data <- filter(data, !(difference == 0))

power_func <- function(x, critical_values){

  redMeth <- unique(x$ReductionMethod)
  ty <- unique(x$type)
  redDim <- unique(x$ReducedDimension)
  pop <- unique(x$populations)
  te <- unique(x$test)
  sa <- unique(x$Samples)

  critical_value <- filter(critical_values, ReductionMethod == redMeth, type == ty,
                           ReducedDimension == redDim, populations == pop, test == te,
                           Samples == sa)$Crit

data.frame(power = mean(x$value >= critical_value))
}

power <- ddply(power_data, .variables = c("ReductionMethod", "type",
                                          "ReducedDimension", "populations",
                                          "test", "difference", "Samples"),
               .fun = power_func, critical_values = critical_values)


