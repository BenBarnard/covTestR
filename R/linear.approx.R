linear.approx <- function (boot.out, L = NULL, index = 1, type = NULL, t0 = NULL,
          t = NULL, ...)
{
  f <- boot.array(boot.out)
  n <- length(f[1, ])
  if ((length(index) > 1L) && (is.null(t0) || is.null(t))) {
    warning("only first element of 'index' used")
    index <- index[1L]
  }
  if (is.null(t0)) {
    t0 <- boot.out$t0[index]
    if (is.null(L))
      L <- empinf(boot.out, index = index, type = type,
                  ...)
  }
  else if (is.null(t) && is.null(L)) {
    warning("input 't0' ignored: neither 't' nor 'L' supplied")
    t0 <- t0[index]
    L <- empinf(boot.out, index = index, type = type, ...)
  }
  else if (is.null(L))
    L <- empinf(boot.out, type = type, t = t, ...)
  tL <- rep(t0, boot.out$R)
  strata <- boot.out$strata
  if (is.null(strata))
    strata <- rep(1, n)
  else strata <- tapply(strata, as.numeric(strata))
  S <- length(table(strata))
  for (s in 1L:S) {
    i.s <- seq_len(n)[strata == s]
    tL <- tL + f[, i.s] %*% L[i.s]/length(i.s)
  }
  as.vector(tL)
}
