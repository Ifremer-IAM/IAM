#' Thanks ThinkR for helping me messing some code
#'
#' @author Maxime Jaunatre
#' @importFrom shiny getShinyOption
#'
#' @noRd
get_golem_options <- function(which = NULL){
  if (is.null(which)){
    getShinyOption("golem_options")
  } else {
    getShinyOption("golem_options")[[which]]
  }
}

#' Thanks ThinkR for helping me messing some code
#'
#' @author Maxime Jaunatre
#' @importFrom shiny runApp
#'
#' @noRd
with_golem_options <- function(app, golem_opts){
  app$appOptions$golem_options <- golem_opts
  return(runApp(app))
}


#' Define UI
#'
#' @author Maxime Jaunatre
#' @import shiny
#' @rawNamespace import(shinyjs, except = c(runExample, alert)) #https://github.com/daattali/shinyjs/issues/127
#' @import shinyWidgets
#'
#' @noRd
app_ui <- function() {
  fluidPage(
    chooseSliderSkin("Flat", color = "#112446"),
    useShinyjs(),
    sidebarPanel(
      radioGroupButtons(
        inputId = "bioeco",
        label = h3("Select variable and filter the dataset"),
        choices = c("Biologic",
                    "Economic"),
        justified = TRUE
      ),
      mod_BioplotUI("bio"),
      mod_EcoplotUI("eco"),
      pickerInput(
        inputId = "sim_name", label = "Scenarii",
        choices = "statu_quo", selected = "statu_quo", multiple = TRUE,
        options = list(`actions-box` = TRUE)
      ),
      hr(),
      checkboxInput("chkribbon", "Display quantile envelopp", value = TRUE),
      sliderInput(
        inputId = "time",
        label = "Time period",
        min = min(get_golem_options("input")$year),
        max = max(get_golem_options("input")$year),
        value = range(get_golem_options("input")$year),
        dragRange = TRUE, step = 1, sep = ""
      ),
      actionButton("done", "Return dataset"),
      actionButton("cancel", "Cancel")
    ),
    mainPanel(
      # h1("Here is the beautiful plot"),
      plotOutput("plot")
    )
  )
}

#' Define server function
#'
#' Depend du petit x reactif qui contient:
#'  bioeco : tri sur le fait d'utiliser le module bio ou eco.
#'  sim_name : tri sur les scenarios a presenter
#'  checkribbon : option de representation des effets stockastiques
#'  time : tri sur la periode de simulation a representer
#'
#'  il contient aussi des variables de retour comme
#'  var : permet d'exporter le tableau trie
#'  bio_plot : permet de representer le plot
#'  eco_plot : pas affecte par cette fonciton
#'
#' @import shiny
#' @import Cairo
#' @importFrom ggplot2 ggplot
#' @rawNamespace import(shinyjs, except = c(runExample, alert)) #https://github.com/daattali/shinyjs/issues/127
#'
#'
#' @author Maxime Jaunatre
#'
#' @noRd
app_server <- function(input, output, session) {
  options(shiny.usecairo=TRUE)

  x <- reactiveValues(
    var = NULL,
    bioeco = NULL, sim_name = NULL, chkribbon = NULL, time = NULL,
    bio_plot = ggplot(NULL), eco_plot = ggplot(NULL)
  )

  observeEvent(input$bioeco, {
    x$bioeco <- input$bioeco
  })
  observeEvent(input$sim_name, {
    x$sim_name <- input$sim_name
  })
  observeEvent(input$chkribbon, {
    x$chkribbon <- input$chkribbon
  })
  observeEvent(input$time, {
    x$time <- input$time
  })

  mod_Bioplot_serv("bio", x)
  mod_Ecoplot_serv("eco", x)

  observe({
    output$plot <- renderPlot({

      return(switch(input$bioeco,
                    Biologic = x$bio_plot,
                    Economic = x$eco_plot))
    })
  })

  # End of app ####
  observeEvent(input$done, {
    # cat("You're really not going to like it\n")
    # returnValue <- x$plot
    returnValue <- x$var
    # TODO : add a print or attribute that explain the filters !
    stopApp(returnValue)
  })
  observeEvent(input$cancel, {
    # cat("Life? Don't talk to me about life.\n")
    stopApp()
  })
}


#' Module shiny pour l'explo graphique.
#'
#' Permet de faire les plots de Fbar, SSB et L_et en fonction du temps.
#'
#' @author Maxime Jaunatre
#' @importFrom shiny shinyApp
#'
#' @noRd
IAM.test_plot <- function(object) {

  if(inherits(object, "iam_formtbl")){
    warning(paste("object must be a table created by IAM.format_quant().",
                  "This function will be called with default probs values."))
    object <- IAM.format_quant(object)
  }
  if(!inherits(object, "iam_quantbl")){
    stop("object must be a table created by IAM.format_quant().")
  }

  Biovars <- c(
    "Fbar", "SSB", "L_et" # "ratio_L_et_tac",
    # (N = courbe d'age)
  )
  Ecovars <- c(
    "nbv_f", "effort2_f", "GVLav_f", "gva_f", "gp_f", "wageg_f",  "wagen_f"
  )

  dinsp <- unique(object$species[!is.na(object$species) & !is.na(object$age)])
  if(length(dinsp) == 0){
    dinsp <- NULL
  }

  res <- with_golem_options(
    app = shinyApp(ui = app_ui, server = app_server),
    golem_opts = list(input = object, Biovars = Biovars, Ecovars = Ecovars,
                      dinsp = dinsp)
  )
  if(is.null(res)){
    return(invisible())
  } else {
    return(res)
  }
}
