bhat_func <- function(ahat21, ahat22){
  ahat21 / ahat22
}

ahat_func <- function(n, p, sample_cov){
  (((n - 1) ^ 2) / (p * (n -2) * (n + 1))) * (tr(sample_cov ^ 2) - (1 / (n - 1)) * (tr(sample_cov)) ^ 2)
}

ahatStar4_func <- funciton(tau, p, S, n){
  (tau / p) *
    (tr(S ^ 4) +
       (-4 / n) * tr(S ^ 3) * tr(S) +
       (-((2 * n ^ 2) + 3 * n - 6) / (n * ((n ^ 2) + n + 2))) * (tr(S ^ 2) ^ 2) +
       ((2 * (5 * n + 6)) / (n * ((n ^ 2) + n + 2))) * tr(S ^ 2) * (tr(S) ^ 2) +
       (-(5 * n + 6) / ((n ^ 2) * ((n ^ 2) + n + 2))) * tr(S) ^ 4)
}

overall_cov_func <- function(A1, A2, n1, n2){
  (1 / (n1 + n2 - 2)) * (A1 + A2)
}

A_func <- function(data){
  sum((data - mean) %*% t((data - mean)))
}

tau_func <- function(n){
  ((n ^ 5) * ((n ^ 2) + n + 2)) / ((n + 1) * (n + 2) * (n + 4) * (n + 6) * (n - 1) * (n - 3))
}

ahat2_func <- function(n, p, S){
  ((n ^ 2) / (p * (n - 1) * (n + 2))) * (tr(S ^ 2) - (1 / n) * (tr(S) ^ 2))
}

deltahat2_func <- function(astar4, p , ahat2, n){
  4 * (((2 * astar4) / (p * ahat2)) * sum(1 / (n - 1)) + sum(1 / ((n - 1) ^ 2)))
}

chaipitak_test <- function(){

}
