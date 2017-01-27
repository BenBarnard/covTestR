#' Test of Equality of Covariances given by Schott 2007
#'
#' Performs 2 and k sample equality of covariance matrix test using Schott 2007
#'
#' @param x data as data.frame, grouped_df, resample or matrix object
#' @param group group or population variable
#' @param ... other options passed to functions
#'
#' @details Consider the test for the equality of covariance matrices with
#' \deqn{\text{H}_0: \boldsymbol{\Sigma}_1 = \boldsymbol{\Sigma}_2 = \ldots = \boldsymbol{\Sigma}_k}
#' and
#' \deqn{\text{H}_{\text{A}}: \boldsymbol{\Sigma}_1 \neq \boldsymbol{\Sigma}_2 \neq \ldots \neq \boldsymbol{\Sigma}_k,}
#' where \eqn{\boldsymbol{\Sigma_i}} are the population covariance matrix parameters with \eqn{n_i\textbf{S}_i \sim W_k\left( \boldsymbol{\Sigma}_i, n_i \right)}
#'
#' \deqn{T_{Sc} = \sum \limits^{k}_{i < j} \frac{ \left( \hat{a}_{2i_{Sc}} + \hat{a}_{2j_{Sc}} - \frac{2}{p}tr \left( S_i S_j \right) \right) ^ 2}{\theta_{Sc}}}
#'
#'
#' @return Test statistic of the hypothesis test
#'
#' @export
#'
#' @references Schott, J. R. (2007). A test for the equality of covariance matrices when
#' the dimension is large relative tothe sample sizes.Computational Statistics & Data Analysis,
#' 51(12):6535–6542.
#'
#' @examples Schott2007_test(iris, group = Species)
#'
Schott2007_test <- function(x, ...) {
  UseMethod("Schott2007_test")
}

#' @export
#' @rdname Schott2007_test
#' @importFrom lazyeval expr_find
#' @importFrom lazyeval lazy_dots
Schott2007_test.data.frame <- function(x, group, ...){
  dataDftoMatrix(data = x,
                 group = expr_find(group),
                 method = expr_find(Schott2007_test.matrix),
                 .dots = lazy_dots(...))
}

#' @export
#' @rdname Schott2007_test
#' @importFrom lazyeval expr_find
#' @importFrom lazyeval lazy_dots
Schott2007_test.grouped_df <- function(x, ...){
  dataDftoMatrix(data = x,
                 group = attributes(x)$vars[[1]],
                 test = expr_find(Schott2007_test.matrix))
}

#' @export
#' @rdname Schott2007_test
#' @importFrom lazyeval expr_find
#' @importFrom lazyeval lazy_dots
Schott2007_test.resample <- function(x, ...){
  x <- as.data.frame(x)
  dataDftoMatrix(data = x,
                 group = attributes(x)$vars[[1]],
                 method = expr_find(Schott2007_test.matrix),
                 .dots = lazy_dots(...))
}

#' @export
#' @rdname Schott2007_test
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#' @importFrom stats cov
#'
Schott2007_test.matrix<- function(...){
  ls <- lazy_dots(...)
  matrix_ls <- lazy_eval(ls[str_detect(names(ls), "x.")])
  names(matrix_ls) <- str_replace(names(matrix_ls), "x.", "")


    ns <- lapply(matrix_ls, function(matrix){
      nrow(matrix)
    })

    p <- lapply(matrix_ls, function(matrix){
      ncol(matrix)
    })

    A_ls <- lapply(matrix_ls, A_func)

    sample_covs <- lapply(matrix_ls, cov)

    overall_cov <- overall_cov_func(A_ls, ns)

    ahat2i <- mapply(ahat2i_func, ns, p, sample_covs, SIMPLIFY = FALSE)
    ahat2 <- ahat2_func(ns, overall_cov, p[[1]])


  Schott2007(ns, p[[1]], ahat2, ahat2i, sample_covs)
}

#' @keywords internal
#' @importFrom utils combn
Schott2007 <- function(ns, p, ahat2, ahat2i, sample_covs){
  comb <- combn(length(ns), 2, simplify = FALSE)
  theta <- 4 * (ahat2 ^ 2) * (Reduce(`+`, lapply(comb, function(x){
    ((1 / ns[[x[1]]]) + (1 / ns[[x[2]]])) ^ 2})) +
      (length(ns) - 1) * (length(ns) - 2) * Reduce(`+`, lapply(ns, function(x){x ^ 2})))

  Reduce(`+`, lapply(comb, function(x){
    ((ahat2i[[x[1]]] + ahat2i[[x[2]]] -
       (2 / p) * tr(sample_covs[[x[1]]] %*% sample_covs[[x[2]]])) ^ 2) / theta
    }))
}
