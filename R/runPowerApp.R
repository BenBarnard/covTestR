#' Run Power Shiny App
#'
#' @export
#' @importFrom shiny runApp
#'
runPowerApp <- function() {
  appDir <- system.file("power_app", package = "EqualCov")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `EqualCov`.", call. = FALSE)
  }

  runApp(appDir, display.mode = "normal", launch.browser = TRUE)
}
