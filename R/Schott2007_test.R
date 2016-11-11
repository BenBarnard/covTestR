#' Test of Equality of Covariances given by Schott 2007
#'
#' @param data tidy data frame
#' @param ... other
#'
#' @return Test Statistic for Schott 2007
#' @export
#'
#' @examples Schott2007_test(mcSamples(c(0,0,0), diag(1, 3), 10, 3), group = population)
#'
Schott2007_test <- function(data, ...) {
  UseMethod("Schott2007_test")
}

#' @export
#'
#' @importFrom lazyeval expr_find
#'
Schott2007_test.data.frame <- function(x, group, ...){
  dataDftoMatrix(data = x,
                 group = expr_find(group),
                 test = expr_find(Schott2007_test.matrix))
}

#' @export
#'
#' @importFrom lazyeval lazy_dots
#' @importFrom lazyeval lazy_eval
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#'
Schott2007_test.matrix<- function(...){
  ls <- lazy_dots(...)
  matrix_ls <- lazy_eval(ls[str_detect(names(ls), "x.")])
  names(matrix_ls) <- str_replace(names(matrix_ls), "x.", "")

<<<<<<< HEAD
    ni <- lapply(matrix_ls, function(matrix){
      nrow(matrix) - 1
=======
    ns <- lapply(matrix_ls, function(matrix){
      nrow(matrix)
>>>>>>> 8250b7d97832b1a55bc7941f5f470838d46915cd
    })

    p <- lapply(matrix_ls, function(matrix){
      ncol(matrix)
    })

<<<<<<< HEAD
    sample_covs <- lapply(matrix_ls, cov)

    scatterMat <- lapply(matrix_ls, A_func)

    overalln <- Reduce(`+`, ni)

    overallScatter <- Reduce(`+`, scatterMat) / overalln


    ahat2 <- (1 / ((overalln - 1) * (overalln + 2) * p[[1]])) * (tr(overallScatter %*% overallScatter) - (1 / overalln) * (tr(overallScatter)) ^ 2)

    npops <- length(matrix_ls)

    its <- combn(seq(npops), 2, simplify = FALSE)


    firstCSum <- Reduce(`+`, lapply(its, function(x){
      (ci_func(ni[[x[1]]], p[[1]]) + ci_func(ni[[x[2]]], p[[1]])) ^ 2
    }))

    secondCSum <- Reduce(`+`, lapply(ni, function(x){
      (ci_func(x, p[[1]])) ^ 2
    }))

    theta2 <- 4 * (ahat2 ^ 2) * (firstCSum + (npops - 1) * (npops - 2) * secondCSum)


  Schott2007_test.default(theta2, sample_covs)
=======
    A_ls <- lapply(matrix_ls, A_func)

    sample_covs <- lapply(matrix_ls, cov)

    overall_cov <- Reduce(`+`, A_ls) / Reduce(`+`, lapply(ns, function(x){x - 1}))

    ahat2i <- mapply(ahat2i_func, ns, p, sample_covs, SIMPLIFY = FALSE)
    ahat2 <- ahat2_func(ns, overall_cov, p[[1]])


  Schott2007_test.default(ns, p[[1]], ahat2, ahat2i, sample_covs)
>>>>>>> 8250b7d97832b1a55bc7941f5f470838d46915cd
}

#' @export
#'
<<<<<<< HEAD
Schott2007_test.default <- function(theta2, sample_covs){

  its <- combn(seq(length(sample_covs)), 2, simplify = FALSE)

  Reduce(`+`, lapply(its, function(x){
    (tr(sample_covs[[x[1]]] - sample_covs[[x[2]]]) ^ 4) / theta2
    }))

=======
Schott2007_test.default <- function(ns, p, ahat2, ahat2i, sample_covs){
  comb <- combn(length(ns), 2, simplify = FALSE)
  theta <- 4 * ahat2 * (Reduce(`+`, lapply(comb, function(x){
    ((1 / ns[[x[1]]]) + (1 / ns[[x[2]]])) ^ 2})) +
      (length(ns) - 1) * (length(ns) - 2) * Reduce(`+`, lapply(ns, function(x){x ^ 2})))

  Reduce(`+`, lapply(comb, function(x){
    ((ahat2i[[x[1]]] + ahat2i[[x[2]]] -
       (2 / p) * tr(sample_covs[[x[1]]] %*% sample_covs[[x[2]]])) ^ 2) / theta
    }))
>>>>>>> 8250b7d97832b1a55bc7941f5f470838d46915cd
}
