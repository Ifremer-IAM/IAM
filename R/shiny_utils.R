#' Accessing app arguments
#'
#' @param which the name of an argument passed to run_app function. char.
#'
#' @importFrom shiny getShinyOption
#'
#' @details Thanks ThinkR for helping me messing some code.
#' @noRd
get_golem_options <- function(which = NULL){
  if (is.null(which)){
    getShinyOption("golem_options")
  } else {
    getShinyOption("golem_options")[[which]]
  }
}


#' Allow to set golem options for an app.
#'
#' @param app a shiny app
#' @param golem_opts list of options. Options need to be named.
#'
#' @importFrom shiny runApp
#'
#' @details Thanks ThinkR for helping me messing some code
#' @noRd
with_golem_options <- function(app, golem_opts){
  app$appOptions$golem_options <- golem_opts
  return(runApp(app))
}
