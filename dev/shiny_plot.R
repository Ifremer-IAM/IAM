library(shiny)
# library(shiny.i18n)
library(shinyWidgets)
# library(shinythemes)
library(shinyjs)
library(tidyverse)

devtools::load_all()
load("dev/data/inputIFR.RData")
load("dev/data/argumIFR.RData")
sim1984 <- IAM::IAM.model(objArgs = IAM_argum_1984, objInput = IAM_input_1984,
                          verbose = FALSE, force_t = NULL)

# Thanks ThinkR for helping me messing some code
get_golem_options <- function(which = NULL){
  if (is.null(which)){
    getShinyOption("golem_options")
  } else {
    getShinyOption("golem_options")[[which]]
  }
}

with_golem_options <- function(app, golem_opts){
  app$appOptions$golem_options <- golem_opts
  return(runApp(app))
}

# Define UI
app_ui <- function() {
  fluidPage(
    chooseSliderSkin("Flat", color = "#112446"),
    useShinyjs(),
    sidebarPanel(
      h3("Select variable and filter the dataset"),
      pickerInput(
        inputId = "var2rep",
        label = "Variables",
        choices = list(
          Bio = get_golem_options("Biovars"),
          Eco = get_golem_options("Ecovars")
        )
      ),
      hidden(
        pickerInput(
          inputId = "species", label = "Species",
          choices = "Fish", selected = "Fish", multiple = TRUE
        ),
        pickerInput(
          inputId = "fleet", label = "Fleet",
          choices = "Boat", selected = "Boat", multiple = TRUE,
          options = list(`actions-box` = TRUE)
        ),
        pickerInput(
          inputId = "metier", label = "Metier",
          selected = "Net", choices = "Net", multiple = TRUE,
          options = list(`actions-box` = TRUE)
        ),
        pickerInput(
          inputId = "age", label = "Age",
          selected = "Young", choices = "Young", multiple = TRUE,
          options = list(`actions-box` = TRUE)
        ),
        checkboxInput("sumage", "Sum selected ages")
      ),
      useSweetAlert(),
      hidden(
        actionBttn(
          inputId = "missing", label = "Information",
          style = "pill", color = "warning"
        )
      ),
      hr(),
      sliderInput(
        inputId = "time",
        label = "Time period",
        min = min(get_golem_options("input")@specific$times),
        max = max(get_golem_options("input")@specific$times),
        value = range(get_golem_options("input")@specific$times),
        dragRange = TRUE, step = 1, sep = ""
      ),
      actionButton("done", "Return dataset"),
      actionButton("cancel", "Cancel")
    ),
    mainPanel(
      h1("Here is the beautiful plot"),
      verbatimTextOutput(outputId = "desc"),
      plotOutput("plot")
    )
  )
}

# Define server function
app_server <- function(input, output, session) {
  output$desc <- renderText(get_golem_options("input")@desc)

  x <- reactiveValues(var = NULL, fullvar = NULL, plot = NULL)
  def <- reactiveValues(specie = NULL, fleet = NULL, metier = NULL, age = NULL)

  observeEvent(input$var2rep, {
    var <- format_var(input$var2rep, get_golem_options("input"))
    # get columns with some values
    c_axes_opt <- FALSE
    if (is.null(var)) {
      axes <- "Variable not computed."
    } else {
      axes <- colnames(var[, -1])[!apply(var[, -1], 2,
                                         function(x) all(is.na(x)))]
      c_axes_opt <- TRUE
    }
    axes <- axes[!axes %in% c("year", "value")]

    def_axes <- c("species", "fleet", "metier", "age")
    for (axe in def_axes) {
      if (axe %in% axes) {
        tmp <- as.character(unique(var[[axe]]))
      } else {
        tmp <- switch(axe,
          species = "Fish",
          fleet = "Boat",
          metier = "Net",
          age = "Young"
        )
      }
      updatePickerInput(session, inputId = axe, choices = tmp, selected = tmp)

      if(axe == "age"){
        updateCheckboxInput(session, "sumage", value = FALSE)
        toggle("sumage", condition = "age" %in% axes)
      }

      def[[axe]] <- tmp
      toggle(axe, condition = axe %in% axes)
    }

    toggle("missing", condition = is.null(var))
    x$fullvar <- var
  })

  observeEvent(input$missing, {# should not trigger.
    sendSweetAlert(
      session = session,
      title = "Information",
      text = paste(
        "This variable was not computed by this simulation or",
        "lack every dimension available usually."
      ),
      type = "Warning"
    )
  })

  observeEvent(input$species, {def$species <- input$species })
  observeEvent(input$fleet, {def$fleet <- input$fleet })
  observeEvent(input$metier, {def$metier <- input$metier })
  observeEvent(input$age, {def$age <- input$age })

  observe({
    N_sp <- !is.null(def$species)
    N_fl <- !is.null(def$fleet)
    N_me <- !is.null(def$metier)
    N_ag <- !is.null(def$age)
    N_sumage <- input$sumage

    # def$species <- switch(any(def$species == "Fish"),NULL)
    if (any(def$species == "Fish")) def$species <- NULL
    if (any(def$fleet == "Boat")) def$fleet <- NULL
    if (any(def$metier == "Net")) def$metier <- NULL
    if (any(def$age == "Young")) def$age <- NULL

    print(input$sumage)

    if (is.null(x$fullvar)) {
      df <- tibble()
    } else {
      df <- x$fullvar %>%
        `if`(N_sp, filter(., species %in% def$species), .) %>%
        `if`(N_fl, filter(., fleet %in% def$fleet), .) %>%
        `if`(N_me, filter(., metier %in% def$metier), .) %>%
        `if`(N_ag, filter(., age %in% def$age), .) %>%
        `if`(N_sumage,
             group_by(., species, fleet, metier, year) %>%
             summarise(., value = sum(value), .groups = "keep") %>%
             ungroup(.), .) %>%
        filter(year >= input$time[1], year <= input$time[2])

    }

    x$var <- df

    if (nrow(df) > 0) {
      style <- aes(linetype = species, color = age)
      style <- style[c(N_sp, N_ag & !input$sumage)]

      p <- ggplot(df, aes(x = year, y = value)) +
        { if (N_fl & !N_me) facet_grid(fleet ~ .) } +
        { if (!N_fl & N_me) facet_grid(. ~ metier) } +
        { if (N_fl & N_me) facet_grid(fleet ~ metier) } +
        geom_line(style) +
        guides(x = guide_axis(angle = 90)) +
        ggtitle(paste(input$var2rep)) +
        NULL
    } else {
      p <- ggplot(NULL)
    }

    x$plot <- p
    output$plot <- renderPlot({
      p
    })

  })

  # End of app ####
  observeEvent(input$done, {
    cat("You're really not going to like it\n")
    returnValue <- 42
    # returnValue <- x$plot
    returnValue <- x$var
    # TODO : add a print or attribute that explain the filters !
    stopApp(returnValue)
  })
  observeEvent(input$cancel, {
    cat("Life? Don't talk to me about life.\n")
    stopApp()
  })
}

run_app <- function(object) {

  Biovars <- c("F","B", "SSB", "L_et", "Ltot", "Fbar", "N")
  Ecovars <- c("nbv_f", "effort1_f", "effort2_f")

  res <- with_golem_options(
    app = shinyApp(ui = app_ui, server = app_server),
    golem_opts = list(input = object, Biovars = Biovars, Ecovars = Ecovars)
  )
  if(is.null(res)){
    return(invisible())
  } else {
    return(res)
  }
}

plot <- run_app(sim1984)


format_var("F", sim1984)

