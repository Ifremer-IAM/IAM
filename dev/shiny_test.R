library(shiny)
library(shiny.i18n)
library(shinyWidgets)
library(shinythemes)
library(shinyjs)


a <- c(
  "A1" = "A1",
  "A2" = "A2",
  "A3" = "A3",
  "A4" = "A4"
)

sp <- c("HKE" = "HKE", "MUT" = "MUT")

ALLVarRep <- c("B","SSB","Ctot","Ytot","Yfmi","Ffmi","Zeit","Fbar","Foth","mu_nbds","mu_nbv","N","Ystat","Lstat","Dstat","Eff",
               "GVL_fme","StatGVL_fme","GVLtot_fm","GVLav_f","vcst_fm","vcst_f","rtbs_f","gp_f","ps_f","gcf_f","gva_f","cs_f","sts_f","rtbsAct_f",
               "csAct_f","gvaAct_f","gcfAct_f","psAct_f","stsAct_f","ccwCr_f","GVLtot_f","wagen_f","L_efmit","D_efmit",
               "Fr_fmi","C_efmit","P","Pstat")

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
        checkboxInput(ns("sr_mod"), "Stock Recrutement Model", value = TRUE),
        selectInput(ns("rec_typ"), "Type",
                    c(
                      "Mean" = "mean",
                      "Hockey-Stick" = "Hockey-Stick",
                      "Beverton-Holt" = "Beverton-Holt",
                      "Ricker" = "Ricker",
                      "Shepherd" = "Shepherd",
                      "Quadratic-HS" = "Quadratic-HS",
                      "Smooth-HS" = "Smooth-HS"
                    ),
                    selectize = FALSE
        ),
        numericInput(ns("a"), label = "a parameter", value = 1),
        numericInput(ns("b"), label = "b parameter", value = 0),
        numericInput(ns("c"), label = "c parameter", value = 0),
        awesomeRadio(
          inputId = ns("Noise_dist"), label = "Noise dist.",
          choices = c("Norm","LogN"), selected = "Norm", inline = TRUE
        ),
        numericInput(ns("noise_dev"), label = "Noise st.dev.", value = 0)
      ),
      column(
        width = 4,
        checkboxInput(ns("stocksim"), "StockSim", value = FALSE),
        awesomeRadio(inputId = ns("sim_typ"), label = "Type",
                     choices = c("1","2","3"), selected = "1"
        )
      )
    )
  )
}

# Define UI
ui <- fluidPage(
  useShinyjs(),
  fluidRow(
    column(
      width = 4,
      div( # Recruitement panels ####
        class = "option-group",
        tabsetPanel(id = "tabs")
      ),
      div( # Iterative panel ####
        class = "option-group",
        checkboxInput("iter", "Iterative", value = FALSE),
        numericInput("niter", label = "Number of iteration", value = 500),
        multiInput(
          inputId = "var_rep",
          label = "Output variables",
          choices = NULL,
          choiceNames = ALLVarRep,
          choiceValues = ALLVarRep
        )
      )
    ),
    column(
      width = 8,
      fluidRow(
        column(
          width = 6,
          div( # Management panel ####
            class = "option-group",
            checkboxInput("manag", "Management", value = FALSE),
            selectInput("man_contr", "Control",
                        c("Nb vessels"="Nb vessels","Nb trips"="Nb trips"),
                        selectize = FALSE
            ),
            selectInput("variables", "Output variables",
                        c("TAC"="TAC","Fbar"="Fbar","TAC->Fbar"="TAC->Fbar"),
                        selectize = FALSE
            ),
            awesomeRadio(
              inputId = "man_typ", label = "Type",choices = c("+","x"),
              selected = "+", inline = TRUE
            ),
            awesomeRadio(
              inputId = "man_upd", label = "Update",choices = c("Yes","No"),
              selected = "Yes", inline = TRUE
            ),
            numericInput("up_bound", label = "Upper bound", value = 0),
            numericInput("low_bound", label = "Lower boun", value = 0)
          )
        ),
        column(
          width = 6,
          div( # Economic panel ####
            class = "option-group",
            checkboxInput("eco", "Economic", value = FALSE),
            # conditionalPanel("input.eco",
            radioGroupButtons(
              inputId = "Eco_typ", choices = c("Complete","DCF"),
              justified = TRUE
            ),
            awesomeRadio(
              inputId = "cc", label = "Adj", choices = c("1","2"),
              selected = "1", inline = TRUE
            ),
            awesomeRadio(
              inputId = "ss", label = "ue_choice", choices = c("1","2"),
              selected = "1", inline = TRUE
            ),
            awesomeRadio(
              inputId = "ii", label = "oths", choices = c("0", "1"),
              selected = "1", inline = TRUE
            ),
            awesomeRadio(
              inputId = "kk", label = "othsFM", choices = c("0", "1"),
              selected = "1", inline = TRUE
            )
            # )
          )
        )
      ),
      div(
        class = "option-group",
        checkboxInput("scenar", "Scenario", value = FALSE),
        selectInput("scen_var", "Output variables",
          a,
          selectize = FALSE
        )
      )
    )
  )
)

# Define server function
server <- function(input, output) {

  for( i in sp){
    appendTab(inputId = "tabs",
              tabPanel(i, mod_spInput(i)) # TODO : need to reload interface with default once tab is changed
    )
  }


  observeEvent(input$iter, {
    toggleState(id = "niter")
    toggleState(id = "var_rep")
  })

  observeEvent(input$eco, {
    toggleState(id = "Eco_typ")
    toggleState(id = "cc")
    toggleState(id = "ss")
    toggleState(id = "ii")
    toggleState(id = "kk")
  })

  observeEvent(input$manag, {
    toggleState(id = "man_contr")
    toggleState(id = "variables")
    toggleState(id = "man_typ")
    toggleState(id = "man_upd")
    toggleState(id = "up_bound")
    toggleState(id = "low_bound")
  })

  observeEvent(input$scenar, {
    toggleState(id = "scen_var")
  })
}



# Create Shiny object
shinyApp(ui = ui, server = server)
