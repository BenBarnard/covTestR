#' Title
#'
#' @param x
#' @param test
#' @param type
#'
#' @return
#' @export
#'
#' @importFrom stringr str_split
#' @importFrom readr read_csv
#' @importFrom readr write_csv
#' @importFrom plyr ddply
#' @importFrom plyr ldply
#' @importFrom plyr l_ply
#' @importFrom stringr str_detect
#'
#'
#' @examples tests(directory = "~/Documents/R/Dissertation/data",
#'                 test_funcs = c("Schott2007_test", "Srivastava2007_test"),
#'                 dimensions = seq(20, 5))
tests <- function(directory, test_funcs, dimensions){

  files <- list.files(directory)

  l_ply(files, function(file, dimensions, test_funcs){
    filegroup <- as.data.frame(str_split(file, " ", simplify = TRUE))
    names(filegroup) <- c("dimension", "samples", "difference", "type", "extension")
    data <- read_csv(paste0(directory, "/", file))

    write_csv(ddply(.data = data, .variables = c("replication", "originaldimensions",
                                       "difference", "ReductionMethod"),
          .fun = function(data, dimensions, test_funcs){
            ldply(.data = dimensions, .fun = function(dimensions, data, test_funcs){
              ldply(lapply(test_funcs, function(x){
                data.frame(replication = unique(data$replication),
                           originaldimensions = unique(data$originaldimensions),
                           difference = unique(data$difference),
                           ReductionMethod = unique(data$ReductionMethod),
                           type = unique(data$type),
                           ReducedDimension = dimensions,
                           populations = max(data$population),
                           test = x,
                           value = do.call(x, list(x = data[, c(seq(dimensions), which(names(data) == "population"))],
                                                   group = quote(population))))
              }))
            }, data = data, test_funcs = test_funcs)
          }, dimensions = c(dimensions, sum(str_detect(names(data), "V"))), test_funcs = test_funcs),
  paste0("~/Documents/R/Dissertation/tests/", file))
  }, dimensions = dimensions, test_funcs = test_funcs)

}
