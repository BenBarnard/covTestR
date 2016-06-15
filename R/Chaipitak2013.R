#' Title
#'
#' @param ahat21
#' @param ahat22
#'
#' @return
#' @export
#'
#' @examples
bhat_func <- function(ahat21, ahat22){
  ahat21 / ahat22
}

#' Title
#'
#' @param n
#' @param p
#' @param sample_cov
#'
#' @return
#' @export
#'
#' @examples
ahat_func <- function(n, p, sample_cov){
  (((n - 1) ^ 2) / (p * (n -2) * (n + 1))) *
    (tr(sample_cov ^ 2) - (1 / (n - 1)) * (tr(sample_cov)) ^ 2)
}

#' Title
#'
#' @param tau
#' @param p
#' @param S
#' @param n
#'
#' @return
#' @export
#'
#' @examples
ahatStar4_func <- function(tau, p, S, n){
  (tau / p) *
    (tr(S ^ 4) +
       (-4 / n) * tr(S ^ 3) * tr(S) +
       (-((2 * n ^ 2) + 3 * n - 6) / (n * ((n ^ 2) + n + 2))) * (tr(S ^ 2) ^ 2) +
       ((2 * (5 * n + 6)) / (n * ((n ^ 2) + n + 2))) * tr(S ^ 2) * (tr(S) ^ 2) +
       (-(5 * n + 6) / ((n ^ 2) * ((n ^ 2) + n + 2))) * tr(S) ^ 4)
  }

#' Title
#'
#' @param A1
#' @param A2
#' @param n1
#' @param n2
#'
#' @return
#' @export
#'
#' @examples
overall_cov_func <- function(A1, A2, n1, n2){
  (1 / (n1 + n2 - 2)) * (A1 + A2)
}

#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
A_func <- function(data){
  sum((data - mean) %*% t((data - mean)))
}

#' Title
#'
#' @param n
#'
#' @return
#' @export
#'
#' @examples
tau_func <- function(n){
  ((n ^ 5) * ((n ^ 2) + n + 2)) / ((n + 1) * (n + 2) * (n + 4) * (n + 6) * (n - 1) * (n - 3))
}

#' Title
#'
#' @param n
#' @param p
#' @param S
#'
#' @return
#' @export
#'
#' @examples
ahat2_func <- function(n, p, S){
  ((n ^ 2) / (p * (n - 1) * (n + 2))) * (tr(S ^ 2) - (1 / n) * (tr(S) ^ 2))
}

#' Title
#'
#' @param astar4
#' @param p
#' @param ahat2
#' @param n
#'
#' @return
#' @export
#'
#' @examples
deltahat2_func <- function(astar4, p , ahat2, n){
  4 * (((2 * astar4) / (p * ahat2)) * sum(1 / (n - 1)) + sum(1 / ((n - 1) ^ 2)))
}

#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr summarise
#' @importFrom dplyr select
#' @importFrom plyr dlpy
#' @importFrom reshape2 acast
#'
#'
#' @examples
chaipitak_test <- function(data){
  browser()
  n_i <- data %>% dlply(c("Group"),
                        function(data){
                          data %>%
                            select(Subjects) %>%
                            unique %>%
                            nrow
                        })

  sample_cov <- data %>% dlply(c("Group"),
                               function(data){
                                 data %>%
                                   select(-Group) %>%
                                   acast(Subjects~Variables,
                                         value.var = "Value")
                               })

  p <- data %>% select(Variables) %>% unique %>% nrow

  A_i <-

}
