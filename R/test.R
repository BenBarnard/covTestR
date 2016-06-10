sim <- replicate(10, simTest(rep(0, 100), diag(1, 100), 100, 3), simplify = TRUE)
sim <- as.data.frame(sim)
ggplot(data = sim) +
  geom_histogram(aes(x = sim, y = ..density..)) +
  stat_function(fun = dnorm, args = list(mean = mean(sqrt(sim$sim))))
 Blah b;ah
