Ei_func <- function(x){
  pn <- pn_func(x)
  sumD2_func(x) / (12 * pn)
}

sumD2_func <- function(x){
n <- nrow(x)
df <- expand.grid(k = seq(1, n),
                  r = seq(1, n),
                  kt = seq(1, n),
                  rt = seq(1, n))
df <- df[!(df$k == df$kt |
             df$k == df$r |
             df$k == df$rt |
             df$r == df$rt |
             df$r == df$kt |
             df$kt == df$rt),]
}

pn_func <- function(x, y = NULL){
  if(y = NULL){
  n <- nrow(x)
  pn <- n * (n - 1) * (n - 2) * (n - 3)
  }else{
    pn <- qn_func(x) * qn_func(y)
  }
  pn
}

qn_func <- function(x){
  n <- nrow(x)
  n * (n - 1)
}

D2ikrktrt_func <- function(x, k, r, kt, rt){
  A2ikrktrt_func(x, k, r, kt, rt) +
    A2ikrktrt_func(x, k, kt, r, rt) +
    A2ikrktrt_func(x, k, rt, kt, r)
}

A2ikrktrt_func <- function(x, k, r, kt, rt){
  t(Dikr_func(x, k, r)) %*% Dikr_func(x, kt, rt) %*%
    t(Dikr_func(x, k, r)) %*% Dikr_func(x, kt, rt)
}

Dikr_func <- function(x, k, r){
  x[k,] - x[r,]
}
