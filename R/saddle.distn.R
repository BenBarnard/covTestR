saddle.distn <- function (A, u = NULL, alpha = NULL, wdist = "m", type = "simp",
          npts = 20, t = NULL, t0 = NULL, init = rep(0.1, d), mu = rep(0.5,
                                                                       n), LR = FALSE, strata = NULL, ...)
{
  call <- match.call()
  if (is.null(alpha))
    alpha <- c(0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2,
               0.5, 0.8, 0.9, 0.95, 0.975, 0.99, 0.995, 0.999)
  if (is.null(t) && is.null(t0))
    stop("one of 't' or 't0' required")
  ep1 <- min(c(alpha, 0.01))/10
  ep2 <- (1 - max(c(alpha, 0.99)))/10
  d <- if (type == "simp")
    1
  else if (is.function(u)) {
    if (is.null(t))
      length(u(t0[1L], ...))
    else length(u(t[1L], ...))
  }
  else 1L + length(u)
  i <- nsads <- 0
  if (!is.null(t))
    npts <- length(t)
  zeta <- matrix(NA, npts, 2L * d - 1L)
  spa <- matrix(NA, npts, 2L)
  pts <- NULL
  if (is.function(A)) {
    n <- nrow(as.matrix(A(t0[1L], ...)))
    if (is.null(u))
      stop("function 'u' missing")
    if (!is.function(u))
      stop("'u' must be a function")
    if (is.null(t)) {
      t1 <- t0[1L] - 2 * t0[2L]
      sad <- saddle(A = A(t1, ...), u = u(t1, ...), wdist = wdist,
                    type = type, d1 = 1, init = init, mu = mu, LR = LR,
                    strata = strata)
      bdu <- bdl <- NULL
      while (is.na(sad$spa[2L]) || (sad$spa[2L] > ep1) ||
             (sad$spa[2L] < ep1/100)) {
        nsads <- nsads + 1
        if (!is.na(sad$spa[2L]) && (sad$spa[2L] > ep1)) {
          i <- i + 1
          zeta[i, ] <- c(sad$zeta.hat, sad$zeta2.hat)
          spa[i, ] <- sad$spa
          pts <- c(pts, t1)
          bdu <- t1
        }
        else bdl <- t1
        if (nsads == npts)
          stop("unable to find range")
        if (is.null(bdl)) {
          t1 <- 2 * t1 - t0[1L]
          sad <- saddle(A = A(t1, ...), u = u(t1, ...),
                        wdist = wdist, type = type, d1 = 1, init = init,
                        mu = mu, LR = LR, strata = strata)
        }
        else if (is.null(bdu)) {
          t1 <- (t0[1L] + bdl)/2
          sad <- saddle(A = A(t1, ...), u = u(t1, ...),
                        wdist = wdist, type = type, d1 = 1, init = init,
                        mu = mu, LR = LR, strata = strata)
        }
        else {
          t1 <- (bdu + bdl)/2
          sad <- saddle(A = A(t1, ...), u = u(t1, ...),
                        wdist = wdist, type = type, d1 = 1, init = init,
                        mu = mu, LR = LR, strata = strata)
        }
      }
      i1 <- i <- i + 1
      nsads <- 0
      zeta[i, ] <- c(sad$zeta.hat, sad$zeta2.hat)
      spa[i, ] <- sad$spa
      pts <- c(pts, t1)
      t2 <- t0[1L] + 2 * t0[2L]
      sad <- saddle(A = A(t2, ...), u = u(t2, ...), wdist = wdist,
                    type = type, d1 = 1, init = init, mu = mu, LR = LR,
                    strata = strata)
      bdu <- bdl <- NULL
      while (is.na(sad$spa[2L]) || (1 - sad$spa[2L] > ep2) ||
             (1 - sad$spa[2L] < ep2/100)) {
        nsads <- nsads + 1
        if (!is.na(sad$spa[2L]) && (1 - sad$spa[2L] >
                                    ep2)) {
          i <- i + 1
          zeta[i, ] <- c(sad$zeta.hat, sad$zeta2.hat)
          spa[i, ] <- sad$spa
          pts <- c(pts, t2)
          bdl <- t2
        }
        else bdu <- t2
        if (nsads == npts)
          stop("unable to find range")
        if (is.null(bdu)) {
          t2 <- 2 * t2 - t0[1L]
          sad <- saddle(A = A(t2, ...), u = u(t2, ...),
                        wdist = wdist, type = type, d1 = 1, init = init,
                        mu = mu, LR = LR, strata = strata)
        }
        else if (is.null(bdl)) {
          t2 <- (t0[1L] + bdu)/2
          sad <- saddle(A = A(t2, ...), u = u(t2, ...),
                        wdist = wdist, type = type, d1 = 1, init = init,
                        mu = mu, LR = LR, strata = strata)
        }
        else {
          t2 <- (bdu + bdl)/2
          sad <- saddle(A = A(t2, ...), u = u(t2, ...),
                        wdist = wdist, type = type, d1 = 1, init = init,
                        mu = mu, LR = LR, strata = strata)
        }
      }
      i <- i + 1
      zeta[i, ] <- c(sad$zeta.hat, sad$zeta2.hat)
      spa[i, ] <- sad$spa
      pts <- c(pts, t2)
      if ((npts%%2) == 0) {
        tt1 <- seq.int(t1, t0[1L], length.out = npts/2 -
                         i1 + 2)[-1L]
        tt2 <- seq.int(t0[1L], t2, length.out = npts/2 +
                         i1 - i + 2)[-1L]
        t <- c(tt1[-length(tt1)], tt2[-length(tt2)])
      }
      else {
        ex <- 1 * (t1 + t2 > 2 * t0[1L])
        ll <- floor(npts/2) + 2
        tt1 <- seq.int(t1, t0[1L], length.out = ll -
                         i1 + 1 - ex)[-1L]
        tt2 <- seq.int(t0[1L], t2, length.out = ll +
                         i1 - i + ex)[-1L]
        t <- c(tt1[-length(tt1)], tt2[-length(tt2)])
      }
    }
    init1 <- init
    for (j in (i + 1):npts) {
      sad <- saddle(A = A(t[j - i], ...), u = u(t[j - i],
                                                ...), wdist = wdist, type = type, d1 = 1, init = init1,
                    mu = mu, LR = LR, strata = strata)
      zeta[j, ] <- c(sad$zeta.hat, sad$zeta2.hat)
      init1 <- sad$zeta.hat
      spa[j, ] <- sad$spa
    }
  }
  else {
    A <- as.matrix(A)
    n <- nrow(A)
    if (is.null(t)) {
      t1 <- t0[1L] - 2 * t0[2L]
      sad <- saddle(A = A, u = c(t1, u), wdist = wdist,
                    type = type, d = d, d1 = 1, init = init, mu = mu,
                    LR = LR, strata = strata)
      bdu <- bdl <- NULL
      while (is.na(sad$spa[2L]) || (sad$spa[2L] > ep1) ||
             (sad$spa[2L] < ep1/100)) {
        if (!is.na(sad$spa[2L]) && (sad$spa[2L] > ep1)) {
          i <- i + 1
          zeta[i, ] <- c(sad$zeta.hat, sad$zeta2.hat)
          spa[i, ] <- sad$spa
          pts <- c(pts, t1)
          bdu <- t1
        }
        else bdl <- t1
        if (i == floor(npts/2))
          stop("unable to find range")
        if (is.null(bdl)) {
          t1 <- 2 * t1 - t0[1L]
          sad <- saddle(A = A, u = c(t1, u), wdist = wdist,
                        type = type, d = d, d1 = 1, init = init,
                        mu = mu, LR = LR, strata = strata)
        }
        else if (is.null(bdu)) {
          t1 <- (t0[1L] + bdl)/2
          sad <- saddle(A = A, u = c(t1, u), wdist = wdist,
                        type = type, d = d, d1 = 1, init = init,
                        mu = mu, LR = LR, strata = strata)
        }
        else {
          t1 <- (bdu + bdl)/2
          sad <- saddle(A = A, u = c(t1, u), wdist = wdist,
                        type = type, d = d, d1 = 1, init = init,
                        mu = mu, LR = LR, strata = strata)
        }
      }
      i1 <- i <- i + 1
      zeta[i, ] <- c(sad$zeta.hat, sad$zeta2.hat)
      spa[i, ] <- sad$spa
      pts <- c(pts, t1)
      t2 <- t0[1L] + 2 * t0[2L]
      sad <- saddle(A = A, u = c(t2, u), wdist = wdist,
                    type = type, d = d, d1 = 1, init = init, mu = mu,
                    LR = LR, strata = strata)
      bdu <- bdl <- NULL
      while (is.na(sad$spa[2L]) || (1 - sad$spa[2L] > ep2) ||
             (1 - sad$spa[2L] < ep2/100)) {
        if (!is.na(sad$spa[2L]) && (1 - sad$spa[2L] >
                                    ep2)) {
          i <- i + 1
          zeta[i, ] <- c(sad$zeta.hat, sad$zeta2.hat)
          spa[i, ] <- sad$spa
          pts <- c(pts, t2)
          bdl <- t2
        }
        else bdu <- t2
        if ((i - i1) == floor(npts/2))
          stop("unable to find range")
        if (is.null(bdu)) {
          t2 <- 2 * t2 - t0[1L]
          sad <- saddle(A = A, u = c(t2, u), wdist = wdist,
                        type = type, d = d, d1 = 1, init = init,
                        mu = mu, LR = LR, strata = strata)
        }
        else if (is.null(bdl)) {
          t2 <- (t0[1L] + bdu)/2
          sad <- saddle(A = A, u = c(t2, u), wdist = wdist,
                        type = type, d = d, d1 = 1, init = init,
                        mu = mu, LR = LR, strata = strata)
        }
        else {
          t2 <- (bdu + bdl)/2
          sad <- saddle(A = A, u = c(t2, u), wdist = wdist,
                        type = type, d = d, d1 = 1, init = init,
                        mu = mu, LR = LR, strata = strata)
        }
      }
      i <- i + 1
      zeta[i, ] <- c(sad$zeta.hat, sad$zeta2.hat)
      spa[i, ] <- sad$spa
      pts <- c(pts, t2)
      if ((npts%%2) == 0) {
        tt1 <- seq.int(t1, t0[1L], length.out = npts/2 -
                         i1 + 2)[-1L]
        tt2 <- seq.int(t0[1L], t2, length.out = npts/2 +
                         i1 - i + 2)[-1L]
        t <- c(tt1[-length(tt1)], tt2[-length(tt2)])
      }
      else {
        ex <- 1 * (t1 + t2 > 2 * t0[1L])
        ll <- floor(npts/2) + 2
        tt1 <- seq.int(t1, t0[1L], length.out = ll -
                         i1 + 1 - ex)[-1L]
        tt2 <- seq.int(t0[1L], t2, length.out = ll +
                         i1 - i + ex)[-1L]
        t <- c(tt1[-length(tt1)], tt2[-length(tt2)])
      }
    }
    init1 <- init
    for (j in (i + 1):npts) {
      sad <- saddle(A = A, u = c(t[j - i], u), wdist = wdist,
                    type = type, d = d, d1 = 1, init = init, mu = mu,
                    LR = LR, strata = strata)
      zeta[j, ] <- c(sad$zeta.hat, sad$zeta2.hat)
      init1 <- sad$zeta.hat
      spa[j, ] <- sad$spa
    }
  }
  pts.in <- (1L:npts)[(abs(zeta[, 1L]) > 1e-06) & (abs(spa[,
                                                           2L] - 0.5) < 0.5 - 1e-10)]
  pts <- c(pts, t)[pts.in]
  zeta <- as.matrix(zeta[pts.in, ])
  spa <- spa[pts.in, ]
  distn <- smooth.spline(qnorm(spa[, 2]), pts)
  quantiles <- predict(distn, qnorm(alpha))$y
  quans <- cbind(alpha, quantiles)
  colnames(quans) <- c("alpha", "quantile")
  inds <- order(pts)
  psa <- cbind(pts[inds], spa[inds, ], zeta[inds, ])
  if (d == 1)
    anames <- "zeta"
  else {
    anames <- rep("", 2 * d - 1)
    for (j in 1L:d) anames[j] <- paste("zeta1.", j, sep = "")
    for (j in (d + 1):(2 * d - 1)) anames[j] <- paste("zeta2.",
                                                      j - d, sep = "")
  }
  dimnames(psa) <- list(NULL, c("t", "gs", "Gs", anames))
  out <- list(quantiles = quans, points = psa, distn = distn,
              call = call, LR = LR)
  class(out) <- "saddle.distn"
  out
}
