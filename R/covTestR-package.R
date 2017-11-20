#' @title Covariance Matrix Testing Functions
#'
#' @description Testing functions for Covariance Matrices. These tests 
#' include high-dimension homogeneity of covariance matrix testing described 
#' by Schott (2007) <doi:10.1016/j.csda.2007.03.004> and high-dimensional 
#' one-sample tests of covariance matrix structure described by 
#' Fisher, et al. (2010) <doi:10.1016/j.jmva.2010.07.004>. Covariance matrix
#' tests use C++ to speed performance and allow larger data sets.
#' @docType package
#' @name covTestR
#' @useDynLib covTestR, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL