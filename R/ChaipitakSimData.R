#' Simulated Data for Chaipitak and Chongaroen 2013
#'
#' @param dim dimensions for covariance matrix
#' @param n number of samples
#' @param covStructure covariance structure
#' @param ... other options passed to functions
#' @param DataMatrix whether to output a data matrix or covariance matrix
#'
#' @return Data either covariance matrix or data matrix
#' @export
#'
#' @importFrom lazyeval lazy_dots
#'
ChaipitakSimData <- function(dim, n, covStructure = "unstructured", ..., DataMatrix = FALSE){
  dots <- lazy_dots(...)
  if(DataMAtrix){
  return()
  }else{
    if("unstructured"){}
    if("compound symmetry"){}
    if("heterogeneous compund symmetry"){}
    if("simple pattern"){}
  }
}
