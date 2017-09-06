#' @title
#' hydroanalyzer_gui
#' @description
#' HydroAnalyzer GUI
#' @import shiny
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
hydroanalyzer_gui <- function() {
  appDir <- system.file("Shiny", "hydroanalyzer", package = "hydroanalyzer")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `hydroanalizer`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
