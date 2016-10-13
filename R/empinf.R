empinf <- function (boot.out = NULL, data = NULL, statistic = NULL, type = NULL,
          stype = NULL, index = 1, t = NULL, strata = rep(1, n), eps = 0.001,
          ...)
{
  if (!is.null(boot.out)) {
    if (boot.out$sim == "parametric")
      stop("influence values cannot be found from a parametric bootstrap")
    data <- boot.out$data
    if (is.null(statistic))
      statistic <- boot.out$statistic
    if (is.null(stype))
      stype <- boot.out$stype
    if (!is.null(boot.out$strata))
      strata <- boot.out$strata
  }
  else {
    if (is.null(data))
      stop("neither 'data' nor bootstrap object specified")
    if (is.null(statistic))
      stop("neither 'statistic' nor bootstrap object specified")
    if (is.null(stype))
      stype <- "w"
  }
  n <- NROW(data)
  if (is.null(type)) {
    if (!is.null(t))
      type <- "reg"
    else if (stype == "w")
      type <- "inf"
    else if (!is.null(boot.out) && (boot.out$sim != "parametric") &&
             (boot.out$sim != "permutation"))
      type <- "reg"
    else type <- "jack"
  }
  if (type == "inf") {
    if (stype != "w")
      stop("'stype' must be \"w\" for type=\"inf\"")
    if (length(index) != 1L) {
      warning("only first element of 'index' used")
      index <- index[1L]
    }
    if (!is.null(t))
      warning("input 't' ignored; type=\"inf\"")
    L <- inf.jack(data, statistic, index, strata, eps, ...)
  }
  else if (type == "reg") {
    if (is.null(boot.out))
      stop("bootstrap object needed for type=\"reg\"")
    if (is.null(t)) {
      if (length(index) != 1L) {
        warning("only first element of 'index' used")
        index <- index[1L]
      }
      t <- boot.out$t[, index]
    }
    L <- empinf.reg(boot.out, t)
  }
  else if (type == "jack") {
    if (!is.null(t))
      warning("input 't' ignored; type=\"jack\"")
    if (length(index) != 1L) {
      warning("only first element of 'index' used")
      index <- index[1L]
    }
    L <- usual.jack(data, statistic, stype, index, strata,
                    ...)
  }
  else if (type == "pos") {
    if (!is.null(t))
      warning("input 't' ignored; type=\"pos\"")
    if (length(index) != 1L) {
      warning("only first element of 'index' used")
      index <- index[1L]
    }
    L <- positive.jack(data, statistic, stype, index, strata,
                       ...)
  }
  L
}
