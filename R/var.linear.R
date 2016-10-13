var.linear <- function (L, strata = NULL)
{
  vL <- 0
  n <- length(L)
  if (is.null(strata))
    strata <- rep(1, n)
  else strata <- tapply(seq_len(n), as.numeric(strata))
  S <- length(table(strata))
  for (s in 1L:S) {
    i.s <- seq_len(n)[strata == s]
    vL <- vL + sum(L[i.s]^2/length(i.s)^2)
  }
  vL
}
