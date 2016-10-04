dim_reduce_singular <- function(data, ...){
  UseMethod("dim_reduce_singular")
}


dim_reduce_singular.data.frame <- function(data, ...){
  do.call(dim_reduce_singular.default,
          data %>%
            dlply(.(Group), dataDftoMatrix))
}


dim_reduce_singular.default <- function(...){
  matrix_ls <- list(...)
browser()



}
