library(plyr)
library(stringr)
library(tidyverse)

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

files <- list.files("~/Box Sync/Dissert Data/100/boxesM")

data <- ldply(files, function(x){
  filegroup <- as.data.frame(str_split(x, " ", simplify = TRUE))
  names(filegroup) <- c("dimension", "samples", "difference", "type", "extension")
  data <- read_csv(paste0("~/Box Sync/Dissert Data/100/boxesM/", x))
  cbind(data, data.frame(Samples = rep(filegroup$samples, nrow(data))))
  })

critical_data <- filter(data, difference == 0)

critical_values <- summarise(group_by(critical_data, ReductionMethod, type,
                                      ReducedDimension, populations, Samples, test),
                             Crit = quantile(value, probs = .95))

power_data <- filter(data, !(difference == 0))

power <- ddply(power_data, .variables = c("ReductionMethod", "type",
                                          "ReducedDimension", "populations",
                                          "test", "difference", "Samples"),
               .fun = power_func, critical_values = critical_values)

write_csv(power, "~/Box Sync/Dissert Data/power/15.csv")

power <- read_csv("~/Box Sync/Dissert Data/power/15.csv")

ggplot(data = filter(power,
                     type == "identity",
                     populations == 3,
                     Samples == 15,
                     !(ReductionMethod == "DataConcatScatterBlock"),
                     difference <= .5)) +
  geom_line(aes(x = difference, y = power, color = test)) +
  facet_grid(ReductionMethod ~ ReducedDimension) +
  theme_bw()

ggplot(data = filter(power,
                     type == "toeplitz",
                     populations == 3,
                     Samples == 15,
                     !(ReductionMethod == "DataConcatScatterBlock"),
                     difference <= 1)) +
  geom_line(aes(x = difference, y = power, color = test)) +
  facet_grid(ReductionMethod ~ ReducedDimension)

ggplot(data = filter(power,
                     type == "identity",
                     populations == 2,
                     Samples == 15,
                     !(ReductionMethod == "DataConcatScatterBlock"),
                     difference <= 1)) +
  geom_line(aes(x = difference, y = power, color = test)) +
  facet_grid(ReductionMethod ~ ReducedDimension)

ggplot(data = filter(power,
                     type == "identity",
                     populations == 3,
                     Samples == 15,
                     !(ReductionMethod == "DataConcatScatterBlock"),
                     difference <= 1)) +
  geom_line(aes(x = difference, y = power, color = test)) +
  facet_grid(ReductionMethod ~ ReducedDimension)

