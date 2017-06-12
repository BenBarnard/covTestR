CompoundSymmetry <- function(p, rho){
  diag(1 - rho, p) + rho * rep(1, p) %*% t(rep(1, p))
}

AutoRegressive <- function(p, rho, kappa){
  sig <- matrix(0, nrow = p, ncol = p)
  for(i in 1:p){
    for(j in 1:p){
      sig[i, j] <- kappa * rho ^ (abs(i - j))
    }
  }
  sig
}

Diagonal <- function(p, min, max){
  diag(runif(p, min, max))
}

Unstructured <- function(){

}
