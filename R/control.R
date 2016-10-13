control <- function (boot.out, L = NULL, distn = NULL, index = 1, t0 = NULL,
          t = NULL, bias.adj = FALSE, alpha = NULL, ...)
{
  if (!is.null(boot.out$call$weights))
    stop("control methods undefined when 'boot.out' has weights")
  if (is.null(alpha))
    alpha <- c(1, 2.5, 5, 10, 20, 50, 80, 90, 95, 97.5, 99)/100
  tL <- dL <- bias <- bias.L <- var.L <- NULL
  k3.L <- q.out <- distn.L <- NULL
  stat <- boot.out$statistic
  data <- boot.out$data
  R <- boot.out$R
  f <- boot.array(boot.out)
  if (bias.adj) {
    if (length(index) > 1L) {
      warning("only first element of 'index' used")
      index <- index[1L]
    }
    f.big <- apply(f, 2L, sum)
    if (boot.out$stype == "i") {
      n <- ncol(f)
      i.big <- rep(seq_len(n), f.big)
      t.big <- stat(data, i.big, ...)[index]
    }
    else if (boot.out$stype == "f")
      t.big <- stat(data, f.big, ...)[index]
    else if (boot.out$stype == "w")
      t.big <- stat(data, f.big/R, ...)[index]
    bias <- mean(boot.out$t[, index]) - t.big
    out <- bias
  }
  else {
    if (is.null(t) || is.null(t0)) {
      if (length(index) > 1L) {
        warning("only first element of 'index' used")
        index <- index[1L]
      }
      if (is.null(L))
        L <- empinf(boot.out, index = index, ...)
      tL <- linear.approx(boot.out, L, index, ...)
      t <- boot.out$t[, index]
      t0 <- boot.out$t0[index]
    }
    else {
      if (is.null(L))
        L <- empinf(boot.out, t = t, ...)
      tL <- linear.approx(boot.out, L, t0 = t0, ...)
    }
    fins <- seq_along(t)[is.finite(t)]
    t <- t[fins]
    tL <- tL[fins]
    R <- length(t)
    dL <- t - tL
    bias.L <- mean(dL)
    strata <- tapply(boot.out$strata, as.numeric(boot.out$strata))
    var.L <- var.linear(L, strata) + 2 * var(tL, dL) + var(dL)
    k3.L <- k3.linear(L, strata) + 3 * cum3(tL, dL) + 3 *
      cum3(dL, tL) + cum3(dL)
    if (is.null(distn)) {
      distn <- saddle.distn((t0 + L)/length(L), alpha = (1L:R)/(R +
                                                                  1), t0 = c(t0, sqrt(var.L)), strata = strata)
      dist.q <- distn$quantiles[, 2]
      distn <- distn$distn
    }
    else dist.q <- predict(distn, x = qnorm((1L:R)/(R + 1)))$y
    distn.L <- sort(dL[order(tL)] + dist.q)
    q.out <- distn.L[(R + 1) * alpha]
    out <- list(L = L, tL = tL, bias = bias.L, var = var.L,
                k3 = k3.L, quantiles = cbind(alpha, q.out), distn = distn)
  }
  out
}
