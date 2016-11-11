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
#'
#' @examples tests(directory = "~/Documents/R/Dissertation/data",
#'                 test_funcs = c("Schott2007_test"),
#'                 dimensions = c(20, 19, 18))
tests <- function(directory, test_funcs, dimensions){

  files <- list.files(directory)

  filegroups <- as.data.frame(str_split(files, " ", simplify = TRUE))
  names(filegroups) <- c("dimension", "samples", "difference", "type", "extension")

  crit_files <- files[filegroups$difference == 0]

  lapply(crit_files, function(file){
    data <- read_csv(paste0(directory, "/", file))
    data <- data[,!(names(data) %in% c("difference"))]

    ddply(.data = data,
          .variables = c("ReductionMethod"),
          .fun = function(data, dimensions, test_funcs){
            if(unique(data$ReductionMethod) == "None"){
              difftrace <- ddply(.data = data,
                                 .variables = c("replication"),
                                 .fun = function(x){
                                   Schott2007_test(x[,!(names(x) %in% c("originaldimensions",
                                                                         "ReductionMethod",
                                                                         "type",
                                                                         "replication"))], population)
                                 })

              lapply(test_funcs, function(funcs){

                test <- ddply(.data = data,
                              .variables = c("replication"),
                              .fun = function(x, funcs){
                                do.call(funcs, list(x[,!(names(x) %in% c("originaldimensions",
                                                                         "ReductionMethod",
                                                                         "type",
                                                                         "replication"))], expr_find(population)))
                              }, funcs = funcs)

               control <- contr(test$V1, difftrace$V1, .95)


              })
            }else{

            }
          }, dimensions = dimensions,
          test_funcs = test_funcs)

  })

  file <- readr::read_csv(paste0(type, "/", x))

  samples <- nrow(file) / (max(file$replication) * max(file$population))
  dimensions <- ncol(file[,-c(1,2)])

  readr::write_csv(plyr::ddply(file, plyr::.(replication),
                               function(file){
                                 do.call(test, list(x = file[,-1],
                                                    group = lazyeval::expr_find(population),
                                                    tidy = FALSE))
                                 }),
                   paste(type, "/", "tests", "/", test, "_", x, sep = ""))
}
