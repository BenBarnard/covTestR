#' Simulate Monte Carlo Samples from a Multivariate Normal Distribution
#'
#' @param meanVec Numeric vector of means for the variables.
#' @param covMat Positive-definite symmetric matrix of the variables.
#' @param maxN Numeric scalar of the maximum of the sample sizes in the simulation.
#' @param ... Other variables used in mvnorm function in the mass package.
#'
#' @return Matrix of maxN multivariate normal samples.
#'
#' @export
#'
#' @examples
#' mcSamples(c(0,0,0), diag(1, 3), 10)
mcSamples <- function(meanVec, covMat, maxN, nCovs, ...){
    dplyr::rename(plyr::ldply(stats::setNames(replicate(nCovs,
              stats::setNames(reshape2::melt(MASS::mvrnorm(n = maxN,
                                                           mu = meanVec,
                                                           Sigma = covMat)),
              c('Subjects', 'Variables', 'Value')),
              simplify = FALSE),
              c(1:nCovs))),
              Group = `.id`)
}

#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
schott2007 <- function(data){
  tnmschott2007(data) / sqrt(theta2schott2007(data))
}

#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
theta2schott2007 <- function(data){
    g <- length(unique(data$Group))
    sumneg2 <- sum(dplyr::summarize(dplyr::group_by(data,
                                                    Group),
                                    sumneg2 = (length(unique(Subjects)) ^ -2))$sumneg2)
    4 * (sumlengthmean(data) + (g - 1) * (g - 2) * sumneg2) * (aschott2007(data) ^ 2)
}

#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
sumlengthmean <- function(data){
  combs <- combn(seq(length(unique(data$Group))), 2)
  combs <- as.data.frame(combs[,combs[1,] - combs[2,] == -1])
  df <- data.frame(Group = unique(data$Group), Index = seq(length(unique(data$Group))))
  Reduce(sum, plyr::alply(combs, 2, function(combs, data, df){
    names(combs) <- "Index"
    combs <- dplyr::inner_join(combs, df, by = c("Index" = "Index"))
    n_i <- length(unique(dplyr::filter(data, Group == combs$Group[1])))
    n_j <- length(unique(dplyr::filter(data, Group == combs$Group[2])))
      ((n_i + n_j) / (n_i * n_j)) ^ 2},
    data = data, df = df))
}

#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
aschott2007 <- function(data){
  n <- nSize(data)
  (m ^ (-1)) * ((n + 2) ^ (-1)) * ((n - 1) ^ (-1)) * (tr(Sschott2007(data) ^ 2) - (n ^ (-1)) * (tr(Sschott2007(data)) ^ 2))
}

#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
tnmschott2007 <- function(data){
  combs <- combn(seq(length(unique(data$Group))), 2)
  combs <- as.data.frame(combs[,combs[1,] - combs[2,] == -1])
  df <- data.frame(Group = unique(data$Group), Index = seq(length(unique(data$Group))))
  Reduce(sum, plyr::alply(combs, 2, function(combs, data, df){
    names(combs) <- "Index"
    combs <- dplyr::inner_join(combs, df, by = c("Index" = "Index"))
    datamat <- plyr::dlply(data,
                c("Group"),
                function(data){
                  reshape2::acast(
                    dplyr::select(data, -Group),
                    Subjects~Variables,
                    value.var = "Value")})
    n_i <- length(unique(dplyr::filter(data, Group == combs$Group[1])))
    n_j <- length(unique(dplyr::filter(data, Group == combs$Group[2])))
    neta_i <- (n_i + 2) * (n_i - 1)
    neta_j <- (n_j + 2) * (n_j - 1)
    S_i <- cov(datamat[[combs$Group[1]]])
    S_j <- cov(datamat[[combs$Group[2]]])
    (1 - (n_i - 2) / neta_i) * tr(S_i ^ 2) +
      (1 - (n_j - 2) / neta_j) * tr(S_j ^ 2) -
      2 * tr(S_i * S_j) -
      n_i * (neta_i ^ (-1)) * (tr(S_i) ^ 2) -
      n_j * (neta_j ^ (-1)) * (tr(S_j) ^ 2)
    },data = data, df = df))
}

#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
Sschott2007 <- function(data){
  Reduce(`+`,
         plyr::llply(
           plyr::dlply(data,
                       c("Group"),
                       function(data){
                         reshape2::acast(
                           dplyr::select(data, -Group),
                           Subjects~Variables,
                           value.var = "Value")}),
           function(data){nrow(data) * cov(data)})) / (length(unique(data$Group)) * length(unique(data$Subjects)))
}

#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
nSize <- function(data){
  sum(dplyr::summarise(dplyr::group_by(data, Group), n = length(unique(Subjects)))$n)}

#' Title
#'
#' @param mat
#'
#' @return
#' @export
#'
#' @examples
tr <- function(mat){
  sum(diag(mat))
}


#' Title
#'
#' @param meanVec
#' @param covMat
#' @param maxN
#' @param nCovs
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
simTest <- function(meanVec, covMat, maxN, nCovs, ...){
  data <- mcSamples(meanVec, covMat, maxN, nCovs)
  schott2007(data)
}
