library(shiny)
library(shinyWidgets)
# library(shinythemes)
library(rhandsontable)
library(shinyjs)

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

mod_spInput <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      # radioGroupButtons(
      #   inputId = "recu_mod",
      #   choices = c("Stock Recrutement Model", "StockSim"),
      #   justified = TRUE
      # ),
      column(
        width = 8,
        checkboxInput(
          ns("sr_mod"), "Stock Recrutement Model",
          value = get_golem_options("arg")$Recruitment[[id]]$modSRactive
        ),
        selectInput(
          ns("rec_typ"), "Type",
          c(
            "Mean" = "mean",
            "Hockey-Stick" = "Hockey-Stick",
            "Beverton-Holt" = "Beverton-Holt",
            "Ricker" = "Ricker",
            "Shepherd" = "Shepherd",
            "Quadratic-HS" = "Quadratic-HS",
            "Smooth-HS" = "Smooth-HS"
          ),
          selectize = FALSE,
          selected = get_golem_options("arg")$Recruitment[[id]]$typeMODsr
        ),
        numericInput(
          ns("a"), label = "a parameter",
          value = get_golem_options("arg")$Recruitment[[id]]$parAmodSR),
        numericInput(
          ns("b"), label = "b parameter",
          value = get_golem_options("arg")$Recruitment[[id]]$parBmodSR),
        numericInput(
          ns("c"), label = "c parameter",
          value = get_golem_options("arg")$Recruitment[[id]]$parCmodSR),
        awesomeRadio(
          inputId = ns("Noise_dist"), label = "Noise dist.",
          choices = c("Norm","LogN"), inline = TRUE,
          selected = c("Norm","LogN")[
            get_golem_options("arg")$Recruitment[[id]]$noiseTypeSR
          ]
        ),
        numericInput(
          ns("noise_dev"), label = "Noise st.dev.",
          value = get_golem_options("arg")$Recruitment[[id]]$wnNOISEmodSR)
      ),
      column(
        width = 4,
        checkboxInput(
          ns("stocksim"), "StockSim",
          value = as.logical(
            get_golem_options("arg")$Recruitment[[id]]$simuSTOCHactive
          )
        ),
        awesomeRadio(
          inputId = ns("sim_typ"), label = "Type", choices = c("1","2","3"),
          selected = as.character(
            get_golem_options("arg")$Recruitment[[id]]$typeSIMUstoch
          )
        )
      )
    )
  )
}

mod_spInput_serv <- function(id, recru){

  moduleServer(
    id,
    ## Below is the module function
    function(input, output, session) {

      observeEvent(input$sr_mod, {
        parts <- c("rec_typ", "a", "b", "c", "Noise_dist", "noise_dev")
        for(p in parts){
          toggleState(id = p, condition = input$sr_mod)
        }
        if(input$sr_mod){
          recru$Recruitment[[id]]$typeMODsr <- input$rec_typ
          recru$Recruitment[[id]]$parAmodSR <- input$a
          recru$Recruitment[[id]]$parBmodSR <- input$b
          recru$Recruitment[[id]]$parCmodSR <- input$c
          # Noise
          recru$Recruitment[[id]]$noiseTypeSR <- match(
            input$Noise_dist, c("Norm", "LogN")
          )
          recru$Recruitment[[id]]$wnNOISEmodSR <- input$noise_dev
        }
        recru$Recruitment[[id]]$modSRactive <- as.numeric(input$sr_mod)
      })

      observeEvent(input$rec_typ, {
        recru$Recruitment[[id]]$typeMODsr <- input$rec_typ
      })
      observeEvent(input$a, {
        recru$Recruitment[[id]]$parAmodSR <- input$a
      })
      observeEvent(input$b, {
        recru$Recruitment[[id]]$parBmodSR <- input$b
      })
      observeEvent(input$c, {
        recru$Recruitment[[id]]$parCmodSR <- input$c
      })
      # Noise
      observeEvent(input$rec_typ, {
        recru$Recruitment[[id]]$noiseTypeSR <- match(
          input$Noise_dist, c("Norm", "LogN")
        )
      })
      observeEvent(input$noise_dev, {
        recru$Recruitment[[id]]$wnNOISEmodSR <- input$noise_dev
      })
      # Stock sim
      observeEvent(input$stocksim, {
        toggleState(id = "sim_typ", condition = input$stocksim)
        if(input$stocksim){
          recru$Recruitment[[id]]$typeSIMUstoch <- as.numeric(input$sim_typ)
        }
        recru$Recruitment[[id]]$simuSTOCHactive <- as.numeric(input$stocksim)
      })
      observeEvent(input$sim_typ, {
        recru$Recruitment[[id]]$typeSIMUstoch <- as.numeric(input$sim_typ)
      })
    }
  )
}

# Define UI
app_ui <- function() {
  fluidPage(
    useShinyjs(),
    fluidRow(
      column(
        width = 4,
        div( # Recruitement panels ####
             class = "option-group",
             tabsetPanel(id = "tabs")
        )
      ),
      column(
        width = 8,
        verbatimTextOutput(outputId = "summary"),
        actionButton("done", "Done"),
        actionButton("cancel", "Cancel")
      )
    )
  )
}

# Define server function
app_server <- function(input, output, session) {

  x <- reactiveValues(
    input = get_golem_options("input"),
    Recruitment = get_golem_options("arg")$Recruitment,
  )



  # add species tab.
  sp <- get_golem_options("spe")$Species
  spDyn <- sp[get_golem_options("spe")$Q == 0]
  for( i in spDyn){
    appendTab(
      inputId = "tabs",
      tabPanel(i, mod_spInput(i))
      # TODO : need to reload interface with default once tab is changed
    )
  }
  updateTabsetPanel(session, "tabs", selected = unname(sp[1]))

  # Recruitment # TODO
  # observe({
    # print(names(input))
    for( i in spDyn){
      mod_spInput_serv(i, x)
    }
  # })
  # EoF Recruitment

  output$summary <- renderText({
    s = capture.output(
      summary(x$input)
    )
    s <- paste0(s, collapse = '\n')
    s
  })

  # TODO : rm Debug part
  observe({
    # TODO : only modify input at the very end.
    x$input@arguments$Recruitment <- x$Recruitment
    print(summary(x$input))
  })

  # End of app ####
  observeEvent(input$done, {
    cat("You're really not going to like it\n")
    returnValue <- 42
    # returnValue <- x$input # TODO : return args
    stopApp(returnValue)
  })
  observeEvent(input$cancel, {
    cat("Life? Don't talk to me about life.\n")
    stopApp()
  })
}

run_app <- function(object, AllVarRep) {

  res <- with_golem_options(
    app = shinyApp(ui = app_ui, server = app_server),
    golem_opts = list(input = object, AllVarRep = AllVarRep,
                      spe = object@specific, arg = object@arguments)
  )
  if(is.null(res)){
    return(invisible())
  } else {
    return(res)
  }
}

# just to test it for now
ALLVarRep = c(
  "B", "SSB", "Ctot", "Ytot", "Yfmi", "Ffmi", "Zeit", "Fbar", "Foth", "mu_nbds", "mu_nbv", "N", "Ystat", "Lstat", "Dstat", "Eff",
  "GVL_fme", "StatGVL_fme", "GVLtot_fm", "GVLav_f", "vcst_fm", "vcst_f", "rtbs_f", "gp_f", "ps_f", "gcf_f", "gva_f", "cs_f", "sts_f", "rtbsAct_f",
  "csAct_f", "gvaAct_f", "gcfAct_f", "psAct_f", "stsAct_f", "ccwCr_f", "GVLtot_f", "wagen_f", "L_efmit", "D_efmit",
  "Fr_fmi", "C_efmit", "P", "Pstat"
)

library(IAM)
data("IAM_argum_1984")
run_app(IAM_argum_1984, ALLVarRep)
a <- summary(IAM_argum_1984)

# TODO : modify app so it initialise with argum values, for SR

IAM::IAM.args(IAM_argum_1984)

