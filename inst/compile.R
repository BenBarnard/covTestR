library(plyr)
library(tidyverse)

files <- list.files("~/Box Sync/Dissert Data/100")

data <- ldply(files, function(x){
  read_csv(paste0("~/Box Sync/Dissert Data/100/", x))
})

critical_data <- filter(data, difference == 0)

critical_values <- summarise(group_by(critical_data, ReductionMethod, type,
                                      ReducedDimension, populations, test),
                             Crit = quantile(value, probs = .95))

power_data <- filter(data, !(difference == 0))

power_func <- function(x, critical_values){

  redMeth <- unique(x$ReductionMethod)
  ty <- unique(x$type)
  redDim <- unique(x$ReducedDimension)
  pop <- unique(x$populations)
  te <- unique(x$test)

  critical_value <- filter(critical_values, ReductionMethod == redMeth, type == ty,
                           ReducedDimension == redDim, populations == pop, test == te)$Crit

data.frame(power = mean(x$value >= critical_value))
}

power <- ddply(power_data, .variables = c("ReductionMethod", "type",
                                          "ReducedDimension", "populations",
                                          "test", "difference"),
               .fun = power_func, critical_values = critical_values)


