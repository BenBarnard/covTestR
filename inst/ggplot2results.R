library(ggplot2)
library(dplyr)
library(tidyr)

load("~/Box Sync/Statistics/Srivastava2010Sim/sim1/powertwotest.RData")
resultspower2Srivastava1 <- powertwotest
load("~/Box Sync/Statistics/Srivastava2010Sim/sim1/powerthreetest.RData")
resultspower3Srivastava1 <- powerthreetest
load("~/Box Sync/Statistics/Srivastava2010Sim/powertwotest.RData")
resultspower2SrivastavaSchott <- powertwotest
load("~/Box Sync/Statistics/Srivastava2010Sim/powerthreetest.RData")
resultspower3SrivastavaSchott <- powerthreetest
resultspower2Srivastava <- bind_rows(resultspower2SrivastavaSchott, resultspower2Srivastava1)
resultspower3Srivastava <- bind_rows(resultspower3SrivastavaSchott, resultspower3Srivastava1)
load("~/Box Sync/Statistics/toeplitz/powertwotest.RData")
resultspower2compound <- powertwotest
load("~/Box Sync/Statistics/toeplitz/powerthreetest.RData")
resultspower3compound <- powerthreetest

dim <- 20
ggplot(data = filter(resultspower2Srivastava, dimension == dim)) +
  geom_line(aes(x = SampleSize, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Unstructured Dimension", dim))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Srivastava/", "Unstructured2Dimension", dim, ".pdf"))

dim <- 40
ggplot(data = filter(resultspower2Srivastava, dimension == dim)) +
  geom_line(aes(x = SampleSize, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Unstructured Dimension", dim))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Srivastava/", "Unstructured2Dimension", dim, ".pdf"))

dim <- 60
ggplot(data = filter(resultspower2Srivastava, dimension == dim)) +
  geom_line(aes(x = SampleSize, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Unstructured Dimension", dim))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Srivastava/", "Unstructured2Dimension", dim, ".pdf"))

dim <- 100
ggplot(data = filter(resultspower2Srivastava, dimension == dim)) +
  geom_line(aes(x = SampleSize, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Unstructured Dimension", dim))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Srivastava/", "Unstructured2Dimension", dim, ".pdf"))

dim <- 200
ggplot(data = filter(resultspower2Srivastava, dimension == dim)) +
  geom_line(aes(x = SampleSize, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Unstructured Dimension", dim))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Srivastava/", "Unstructured2Dimension", dim, ".pdf"))

samp <- 10
ggplot(data = filter(resultspower2Srivastava, SampleSize == samp)) +
  geom_line(aes(x = dimension, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Unstructured Sample", samp))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Srivastava/", "Unstructured2Sample", samp, ".pdf"))

samp <- 20
ggplot(data = filter(resultspower2Srivastava, SampleSize == samp)) +
  geom_line(aes(x = dimension, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Unstructured Sample", samp))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Srivastava/", "Unstructured2Sample", samp, ".pdf"))

samp <- 40
ggplot(data = filter(resultspower2Srivastava, SampleSize == samp)) +
  geom_line(aes(x = dimension, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Unstructured Sample", samp))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Srivastava/", "Unstructured2Sample", samp, ".pdf"))

samp <- 60
ggplot(data = filter(resultspower2Srivastava, SampleSize == samp)) +
  geom_line(aes(x = dimension, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Unstructured Sample", samp))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Srivastava/", "Unstructured2Sample", samp, ".pdf"))

dim <- 20
ggplot(data = filter(resultspower3Srivastava, dimension == dim)) +
  geom_line(aes(x = SampleSize, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Unstructured Dimension", dim))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Srivastava/", "Unstructured3Dimension", dim, ".pdf"))

dim <- 40
ggplot(data = filter(resultspower3Srivastava, dimension == dim)) +
  geom_line(aes(x = SampleSize, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Unstructured Dimension", dim))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Srivastava/", "Unstructured3Dimension", dim, ".pdf"))

dim <- 60
ggplot(data = filter(resultspower3Srivastava, dimension == dim)) +
  geom_line(aes(x = SampleSize, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Unstructured Dimension", dim))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Srivastava/", "Unstructured3Dimension", dim, ".pdf"))

dim <- 100
ggplot(data = filter(resultspower3Srivastava, dimension == dim)) +
  geom_line(aes(x = SampleSize, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Unstructured Dimension", dim))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Srivastava/", "Unstructured3Dimension", dim, ".pdf"))

dim <- 200
ggplot(data = filter(resultspower3Srivastava, dimension == dim)) +
  geom_line(aes(x = SampleSize, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Unstructured Dimension", dim))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Srivastava/", "Unstructured3Dimension", dim, ".pdf"))

samp <- 10
ggplot(data = filter(resultspower3Srivastava, SampleSize == samp)) +
  geom_line(aes(x = dimension, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Unstructured Sample", samp))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Srivastava/", "Unstructured3Sample", samp, ".pdf"))

samp <- 20
ggplot(data = filter(resultspower3Srivastava, SampleSize == samp)) +
  geom_line(aes(x = dimension, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Unstructured Sample", samp))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Srivastava/", "Unstructured3Sample", samp, ".pdf"))

samp <- 40
ggplot(data = filter(resultspower3Srivastava, SampleSize == samp)) +
  geom_line(aes(x = dimension, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Unstructured Sample", samp))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Srivastava/", "Unstructured3Sample", samp, ".pdf"))

samp <- 60
ggplot(data = filter(resultspower3Srivastava, SampleSize == samp)) +
  geom_line(aes(x = dimension, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Unstructured Sample", samp))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Srivastava/", "Unstructured3Sample", samp, ".pdf"))









dim <- 20
ggplot(data = filter(resultspower2compound, dimension == dim)) +
  geom_line(aes(x = SampleSize, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Compound Dimension", dim))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Chapitak/", "Compound2Dimension", dim, ".pdf"))

dim <- 40
ggplot(data = filter(resultspower2compound, dimension == dim)) +
  geom_line(aes(x = SampleSize, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Compound Dimension", dim))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Chapitak/", "Compound2Dimension", dim, ".pdf"))

dim <- 60
ggplot(data = filter(resultspower2compound, dimension == dim)) +
  geom_line(aes(x = SampleSize, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Compound Dimension", dim))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Chapitak/", "Compound2Dimension", dim, ".pdf"))

dim <- 100
ggplot(data = filter(resultspower2compound, dimension == dim)) +
  geom_line(aes(x = SampleSize, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Compound Dimension", dim))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Chapitak/", "Compound2Dimension", dim, ".pdf"))

dim <- 200
ggplot(data = filter(resultspower2compound, dimension == dim)) +
  geom_line(aes(x = SampleSize, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Compound Dimension", dim))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Chapitak/", "Compound2Dimension", dim, ".pdf"))

samp <- 10
ggplot(data = filter(resultspower2compound, SampleSize == samp)) +
  geom_line(aes(x = dimension, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Compound Sample", samp))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Chapitak/", "Compound2Sample", samp, ".pdf"))

samp <- 20
ggplot(data = filter(resultspower2compound, SampleSize == samp)) +
  geom_line(aes(x = dimension, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Compound Sample", samp))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Chapitak/", "Compound2Sample", samp, ".pdf"))

samp <- 40
ggplot(data = filter(resultspower2compound, SampleSize == samp)) +
  geom_line(aes(x = dimension, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Compound Sample", samp))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Chapitak/", "Compound2Sample", samp, ".pdf"))

samp <- 60
ggplot(data = filter(resultspower2compound, SampleSize == samp)) +
  geom_line(aes(x = dimension, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Compound Sample", samp))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Chapitak/", "Compound2Sample", samp, ".pdf"))

dim <- 20
ggplot(data = filter(resultspower3compound, dimension == dim)) +
  geom_line(aes(x = SampleSize, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Compound Dimension", dim))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Chapitak/", "Compound3Dimension", dim, ".pdf"))

dim <- 40
ggplot(data = filter(resultspower3compound, dimension == dim)) +
  geom_line(aes(x = SampleSize, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Compound Dimension", dim))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Chapitak/", "Compound3Dimension", dim, ".pdf"))

dim <- 60
ggplot(data = filter(resultspower3compound, dimension == dim)) +
  geom_line(aes(x = SampleSize, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Compound Dimension", dim))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Chapitak/", "Compound3Dimension", dim, ".pdf"))

dim <- 100
ggplot(data = filter(resultspower3compound, dimension == dim)) +
  geom_line(aes(x = SampleSize, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Compound Dimension", dim))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Chapitak/", "Compound3Dimension", dim, ".pdf"))

dim <- 200
ggplot(data = filter(resultspower3compound, dimension == dim)) +
  geom_line(aes(x = SampleSize, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Compound Dimension", dim))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Chapitak/", "Compound3Dimension", dim, ".pdf"))

samp <- 10
ggplot(data = filter(resultspower3compound, SampleSize == samp)) +
  geom_line(aes(x = dimension, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Compound Sample", samp))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Chapitak/", "Compound3Sample", samp, ".pdf"))

samp <- 20
ggplot(data = filter(resultspower3compound, SampleSize == samp)) +
  geom_line(aes(x = dimension, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Compound Sample", samp))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Chapitak/", "CompoundSample", samp, ".pdf"))

samp <- 40
ggplot(data = filter(resultspower3compound, SampleSize == samp)) +
  geom_line(aes(x = dimension, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Compound Sample", samp))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Chapitak/", "Compound3Sample", samp, ".pdf"))

samp <- 60
ggplot(data = filter(resultspower3compound, SampleSize == samp)) +
  geom_line(aes(x = dimension, y = Power, color = Test, group = Test)) +
  theme_bw() +
  ggtitle(paste("Compound Sample", samp))
ggsave(paste0("/Users/ben_barnard/Desktop/power/Chapitak/", "Compound3Sample", samp, ".pdf"))
