solve_frobenius <- function(d2, covA){
  m <- ncol(covA)
  covA + (sqrt(d2 / (m * m)))
}
