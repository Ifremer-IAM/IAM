library(shiny)
library(shiny.i18n)
library(shinyWidgets)
library(shinythemes)
library(shinyjs)

devtools::load_all()
load("dev/data/inputIFR.RData")
load("dev/data/argumIFR.RData")
sim1984 <- IAM::IAM.model(objArgs = argum1984, objInput = input1984,
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
        )
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

  x <- reactiveValues(var = NULL, fullvar = NULL)
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
        tmp <- unique(var[[axe]])
      } else {
        tmp <- switch(axe,
          species = "Fish",
          fleet = "Boat",
          metier = "Net",
          age = "Young"
        )
      }
      updatePickerInput(session, inputId = axe, choices = tmp, selected = tmp)
      def[[axe]] <- tmp
      toggle(axe, condition = axe %in% axes)
    }


    # if("species" %in% axes){
    #   def$species <- unique(var$species)
    #   updatePickerInput(
    #     session, inputId = "species",
    #     choices = unique(var$species), selected = unique(var$species))
    #   show("species")
    # } else {
    #   # def$species <- NULL
    #   updatePickerInput(
    #     session, inputId = "species",
    #     choices = "Fish", selected = "Fish"
    #   )
    #   hide("species")
    # } # module this part!
    # if("fleet" %in% axes){
    #   def$fleet <- unique(var$fleet)
    #   updatePickerInput(
    #     session, inputId = "fleet",
    #     choices = unique(var$fleet), selected = unique(var$fleet))
    #   show("fleet")
    # } else {
    #   def$fleet <- NULL
    #   updatePickerInput(
    #     session, inputId = "fleet",
    #     choices = "Boat", selected = "Boat")
    #   hide("fleet")
    # }
    # if("metier" %in% axes){
    #   def$metier <- unique(var$metier)
    #   updatePickerInput(
    #     session, inputId = "metier",
    #     choices = unique(var$metier), selected = unique(var$metier))
    #   show("metier")
    # } else {
    #   def$metier <- NULL
    #   updatePickerInput(
    #     session, inputId = "metier",
    #     choices = "Net", selected = "Net")
    #   hide("metier")
    # }
    # if("age" %in% axes){
    #   def$age <- unique(var$age)
    #   updatePickerInput(
    #     session, inputId = "age",
    #     choices = unique(var$age), selected = unique(var$age))
    #   show("age")
    # } else {
    #   def$age <- NULL
    #   updatePickerInput(
    #     session, inputId = "age",
    #     choices = "Young", selected = "Young")
    #   hide("age")
    # }

    toggle("missing", condition = is.null(var))
    x$fullvar <- var
  })

  observeEvent(input$missing, {
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

    # def$species <- switch(any(def$species == "Fish"),NULL)
    if (any(def$species == "Fish")) def$species <- NULL
    if (any(def$fleet == "Boat")) def$fleet <- NULL
    if (any(def$metier == "Net")) def$metier <- NULL
    if (any(def$age == "Young")) def$age <- NULL

    if (is.null(x$fullvar)) {
      df <- tibble()
    } else {
      df <- x$fullvar %>%
        `if`(N_sp, filter(., species %in% def$species), .) %>%
        `if`(N_fl, filter(., fleet %in% def$fleet), .) %>%
        `if`(N_me, filter(., metier %in% def$metier), .) %>%
        `if`(N_ag, filter(., age %in% def$age), .) %>%
        filter(year >= input$time[1], year <= input$time[2])
    }

    if (nrow(df) > 0) {
      style <- aes(linetype = species, color = age)
      style <- style[c(N_sp, N_ag)]

      p <- ggplot(df, aes(x = year, y = value)) +
        { if (N_fl & !N_me) facet_grid(fleet ~ .) } +
        { if (!N_fl & N_me) facet_grid(. ~ metier) } +
        { if (N_fl & N_me) facet_grid(fleet ~ metier) } +
        geom_line(style)
    } else {
      p <- ggplot(NULL)
    }

    output$plot <- renderPlot({
      p
    })
    x$var <- df
  })

  # End of app ####
  observeEvent(input$done, {
    cat("You're really not going to like it\n")
    returnValue <- 42
    # returnValue <- x$var
    stopApp(returnValue)
  })
  observeEvent(input$cancel, {
    cat("Life? Don't talk to me about life.\n")
    stopApp()
  })
}

run_app <- function(object) {

  # vars <- c("B","SSB","Ctot","Ytot","Yfmi","Ffmi","Zeit","Fbar","Foth",
  #           "mu_nbds","mu_nbv","N","Eff", "GVL_fme","GVLtot_fm","GVLav_f",
  #           "rtbs_f","gp_f","ps_f","gcf_f","gva_f","cs_f","sts_f","rtbsAct_f",
  #           "csAct_f","gvaAct_f","gcfAct_f","psAct_f","stsAct_f","ccwCr_f",
  #           "GVLtot_f","wagen_f","L_efmit","D_efmit","Fr_fmi","C_efmit",
  #           "vcst_f","vcst_fm","P", "Ystat","Lstat","Dstat","Pstat",
  #           "StatGVL_fme")

  Biovars <- c("SSB", "Fbar", "N", "Eff")
  Ecovars <- c("GVLtot_f_m", "gvaAct_f")

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

run_app(sim1984)


