#' 'iamArgs' objects creator by code
#'
#' @param object  An \code{\link[IAM]{iamInput-class}} or \code{\link[IAM]{iamArgs-class}} object.
#' @param ... desc argument
#'
#' @name IAM.input2args
#'
#' @export
setGeneric("IAM.input2args", function(object,  ...){
  standardGeneric("IAM.input2args")
}
)

#' @param desc Object descriptor (default value : \code{as.character(NA)}).
#' If not provided, copied the description slot of object.
#'
#' @rdname IAM.input2args
#' @export
setMethod("IAM.input2args", signature("iamInput"),function(object, desc= NULL){


  ALLVarRep = c(
    "B", "SSB", "Ctot", "Ytot", "Yfmi", "Ffmi", "Zeit", "Fbar", "Foth",
    "mu_nbds", "mu_nbv", "N", "Ystat", "Lstat", "Dstat", "Eff",
    "GVL_fme", "StatGVL_fme", "GVLtot_fm", "GVLav_f", "vcst_fm", "vcst_f",
    "rtbs_f", "gp_f", "ps_f", "gcf_f", "gva_f", "cs_f", "sts_f", "rtbsAct_f",
    "csAct_f", "gvaAct_f", "gcfAct_f", "psAct_f", "stsAct_f", "ccwCr_f",
    "GVLtot_f", "wagen_f", "L_efmit", "D_efmit",
    "Fr_fmi", "C_efmit", "P", "Pstat"
  )

  if(is.null(desc)){ desc <- object@desc }
  # init the arg object with shiny default
  ## Create argum ####
  # TODO : this is where to modify the default for the GUI now !
  init_recru <- function(name, object){
    list(modSRactive = 1,
         typeMODsr = "Mean",
         parAmodSR = unname(object@input[[name]]$N_it0[1]),
         parBmodSR = 0, parCmodSR = 0, wnNOISEmodSR = 0, noiseTypeSR = 1,
         simuSTOCHactive = 0, typeSIMUstoch = 1
    )
  }

  sp <- as.list(object@specific$Species) ; names(sp) <- sp
  spDyn <- sp[object@specific$Q == 0]
  Recruitement <- lapply(spDyn, function(x) init_recru(name = x, object = object))
  rm(sp)

  Replicates <- list(active =0, nbIter =500,
                     SELECTvar = ALLVarRep[c(1:20, 22:38, 43:44)] )

  Scenario <- list(active = 0, ALLscenario = names(object@scenario), SELECTscen = 1)

  Eco <- list(
    # active = 0, type = 1,
    # adj = 1, ue_choice = 1, oths = 0, othsFM = 0, # useless
    perscCalc = 0,
    # report = 0, # useless
    dr = 0
  )

  init_gest <- function(object){
    years <- object@specific$times
    fleets <- object@specific$Fleet
    spp <- c(na.omit(object@specific$Species))
    sppStat <- c(na.omit(object@specific$StaticSpp))
    Espece <- c(na.omit(c(spp,sppStat)))[1]

    tacfbar <- matrix(as.numeric(NA),nrow=2,ncol=length(years),
                      dimnames=list(c("TAC","Fbar"),years))
    mfm <- object@input$Fleet$nbv_f_m ; mfm[!is.na(mfm)] <- 1

    Gestion <- list(
      active = 0, control = "Nb vessels", target = "TAC",
      espece = Espece, delay = 2, typeG = 0, upd = 1,
      sup = 0, inf =0,
      tac = tacfbar["TAC",],
      fbar = tacfbar["Fbar",],
      # othSpSup = tabO,
      effSup = matrix(as.numeric(NA),nrow=length(fleets),
                      ncol=length(years),dimnames=list(fleets,years)),
      mfm = mfm
    )

    attr(Gestion$tac, "DimCst") <- as.integer(c(0,0,0,length(Gestion$tac)))
    attr(Gestion$fbar, "DimCst") <- as.integer(c(0,0,0,length(Gestion$fbar)))
    attr(Gestion$effSup, "DimCst") <- as.integer(c(nrow(Gestion$effSup),0,0,
                                                   ncol(Gestion$effSup)))
    attr(Gestion$mfm, "DimCst") <- as.integer(c(dim(Gestion$mfm),0,0))
    return(Gestion)
  }
  Gestion <- init_gest(object)

  argum <- list(Recruitment = Recruitement, Replicates = Replicates,
                Scenario = Scenario, Gestion = Gestion, Eco = Eco)
  ## Copy specific ####
  specif <- object@specific
  args <- new("iamArgs", desc = desc, arguments = argum, specific = specif)

  return(args)

})


#' Adding tabset for each species
#'
#' @import shiny
#' @rawNamespace import(shinyjs, except = c(runExample, alert)) #https://github.com/daattali/shinyjs/issues/127
#' @import shinyWidgets
#'
#' @param id names of the species and so name of the tab
#'
#' @author Maxime Jaunatre
#' @noRd
#'
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

#' @noRd
#' @param recru reactive value to access recruitment slot and modify it.
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

#' @title Shiny app for \code{IAM.args()}.
#'
#' @author Maxime Jaunatre
#'
#' @import shiny
#' @rawNamespace import(shinyjs, except = c(runExample, alert)) #https://github.com/daattali/shinyjs/issues/127
#' @import shinyWidgets
#' @import rhandsontable
#'
#' @noRd
args_app_ui <- function() {
  fluidPage(
    tags$head(
      tags$style(HTML("
      .option_group {
          border: 1px double gray;
          border-radius: 10px !important;
          padding: 5px;
      }
    "))),

    useShinyjs(),
    fluidRow(
      column(
        style='padding:5px;',
        width = 4,
        tags$div(
          class = "option_group", # Recruitement panels ####
          h3("Recruitment"),
          tabsetPanel(id = "tabs")
        )
      ),
      column(
        style='padding:5px;',
        width = 4,
        tags$div(
          class = "option_group", # Management panel ####
          fluidRow(
            column(6, h3("Management")),
            column(6, checkboxInput(
              "manag", "activate",
              value = as.logical(get_golem_options("arg")$Gestion$active)
            ))
          ),
          fluidRow(
            column(6, selectInput(
              "espece", "Species",
              get_golem_options("spe")$Species, # TODO : only dynamic or all species ?
              selectize = FALSE,
              selected =get_golem_options("arg")$Gestion$espece
            )),
            column(6, selectInput(
              "control", "Control",
              c("Nb vessels"="Nb vessels","Nb trips"="Nb trips"),
              selectize = FALSE,
              selected = get_golem_options("arg")$Gestion$control
            ))
          ),
          h4("Target"),
          fluidRow(
            column(5, selectInput(
              "target", label = NULL, #"Target",
              c("TAC"="TAC","Fbar"="Fbar","TAC->Fbar"="TAC->Fbar"),
              selectize = FALSE,
              selected = get_golem_options("arg")$Gestion$target
            )),
            column(7, dropdown(
              inputId = "tacfbar",
              rHandsontableOutput("tac", width = 450),
              width = "500px", up = FALSE, right = TRUE,
              label = "TAC/Fbar", tooltip = TRUE
            ))
          ),
          fluidRow(
            column(6, awesomeRadio(
              inputId = "typeG", label = "Type",choices = c("+","x"),
              inline = TRUE,
              selected = c("+", "x")[get_golem_options("arg")$Gestion$typeG +1]
            )),
            column(6, awesomeRadio(
              inputId = "upd", label = "Update",choices = c("Yes","No"),
              inline = TRUE,
              selected = c("Yes", "No")[get_golem_options("arg")$Gestion$upd]
            ))
          ),
          sliderInput(
            "delay", label = "Delay", step = 1, min = 1,
            max = get_golem_options("spe")$NbSteps,
            value = get_golem_options("arg")$Gestion$delay
          ),
          fluidRow(
            column(6, numericInput(
              "sup", label = "Upper bound",
              value = get_golem_options("arg")$Gestion$sup
            )),
            column(6, numericInput(
              "inf", label = "Lower bound",
              value = get_golem_options("arg")$Gestion$inf
            ))
          ),
          fluidRow(
            column(6, dropdown(
              inputId = "effsup_mat",
              rHandsontableOutput("eff", width = 450),
              width = "500px", up = TRUE,
              label = "eff table", tooltip = FALSE
            )),
            column(6, dropdown(
              inputId = "mfm_mat",
              rHandsontableOutput("mfm", width = 450),
              width = "500px",up = TRUE, right = TRUE,
              label = "mfm table", tooltip = FALSE
            )),
          )
        ),
        tags$div(
          class = "option_group", # Scenar ####
          fluidRow(
            column(6, h3("Scenario")),
            column(6, checkboxInput(
              "scenar", "activate",
              value = as.logical(get_golem_options("arg")$Scenario$active)
            ))
          ),
          selectInput(
            "scen_var", "Output variables",
            get_golem_options("arg")$Scenario$ALLscenario,
            selectize = FALSE,
            selected =  get_golem_options("arg")$Scenario$ALLscenario[
              get_golem_options("arg")$Scenario$SELECTscen
            ]
          )
        )
      ),
      column(
        style='padding:5px;',
        width = 4,
        tags$div(
          class = "option_group", # Economic panel ####
          h3("Economic"),
          awesomeRadio(
            inputId = "perscCalc", label = "perscCalc",
            choices = c("0", "1", "2", "3", "4"), inline = TRUE,
            selected = as.character(get_golem_options("arg")$Eco$perscCalc)
          ),
          numericInput(
            "dr", label = "Discound rate", step = 0.01,
            value = get_golem_options("arg")$Eco$dr
          )
        ),
        tags$div(
          class = "option_group", # Replicates panel ####
          fluidRow(
            column(6, h3("Replicates")),
            column(6, checkboxInput(
              "iter", "activate",
              value = as.logical(get_golem_options("arg")$Replicates$active)
            ))
          ),
          numericInput(
            "nbIter", label = "Number of iteration",
            value = get_golem_options("arg")$Replicates$nbIter
          ),
          multiInput(
            inputId = "var_rep",
            label = "Output variables",
            choices = NULL,
            choiceNames = get_golem_options("AllVarRep"),
            choiceValues = get_golem_options("AllVarRep")
          )
        )
      )
    ),
    fluidRow(
      column(
        4, offset = 4,
        actionButton("done", "OK"),
        actionButton("cancel", "Cancel")
      )
    )
  )
}

#' @noRd
#' @param input basic shiny argument for server function
#' @param output basic shiny argument for server function
#' @param session basic shiny argument for server function
args_app_server <- function(input, output, session) {

  x <- reactiveValues(
    input = get_golem_options("input"),
    Recruitment = get_golem_options("arg")$Recruitment,
    Replicates = get_golem_options("arg")$Replicates,
    Scenario = get_golem_options("arg")$Scenario,
    Gestion = get_golem_options("arg")$Gestion,
    Eco = get_golem_options("arg")$Eco
  )

  # Editable table with rhandsontable ####
  Gestion_tab <- reactiveValues()
  tac <- reactive({
    if (!is.null(input$tac)) {
      DF = hot_to_r(input$tac)
      rownames(DF) <- rownames(get_golem_options("tac"))
    } else {
      if (is.null(Gestion_tab$TAC)){
        DF = get_golem_options("tac")
      }else{
        DF = Gestion_tab$TAC
      }
    }
    Gestion_tab$TAC = DF
    DF
  })
  maxeff <- reactive({
    if (!is.null(input$eff)) {
      DF = hot_to_r(input$eff)
      rownames(DF) <- rownames(get_golem_options("maxeff"))
    } else {
      if (is.null(Gestion_tab$eff))
        DF = get_golem_options("maxeff")
      else
        DF = Gestion_tab$eff
    }
    Gestion_tab$eff = DF
    DF
  })
  mfm <- reactive({
    if (!is.null(input$mfm)) {
      DF = hot_to_r(input$mfm)
      rownames(DF) <- rownames(get_golem_options("mfm"))
    } else {
      if (is.null(Gestion_tab$mfm))
        DF = get_golem_options("mfm")
      else
        DF = Gestion_tab$mfm
    }
    Gestion_tab$mfm = DF
    DF
  })
  # eof tables

  # add species tab.
  sp <- get_golem_options("spe")$Species
  spDyn <- sp[get_golem_options("spe")$Q == 0]
  for( i in spDyn){
    appendTab(inputId = "tabs",
              tabPanel(i, mod_spInput(i))
    )
  }
  updateTabsetPanel(session, "tabs", selected = unname(sp[1]))

  ## Recruitment ####
  for( i in spDyn){
    mod_spInput_serv(i, x)
  }
  # EoF Recruitment

  ## Replicates ####
  observeEvent(input$iter, {
    toggleState(id = "nbIter", condition = input$iter)
    toggleState(id = "var_rep", condition = input$iter)
    x$Replicates$active <- as.numeric(input$iter)
  })
  observeEvent(input$nbIter, {
    x$Replicates$nbIter <- input$nbIter
  })
  observeEvent(input$var_rep, {
    x$Replicates$SELECTvar <- input$var_rep
  })
  # EoF Replicates

  ## Gestion ####
  observeEvent(input$manag, {
    parts <- c("espece", "control", "target", "typeG", "delay", "upd", "sup",
               "inf", "tacfbar", "effsup_mat", "mfm_mat")
    for(p in parts){
      toggleState(id = p, condition = input$manag)
    }
    x$Gestion$active <- as.numeric(input$manag)
  })
  observeEvent(input$control,{
    x$Gestion$control <- input$control
  })
  observeEvent(input$target,{
    x$Gestion$target <- input$target
  })
  observeEvent(input$espece,{
    x$Gestion$espece <- input$espece
  })
  observeEvent(input$delay,{
    x$Gestion$delay <- input$delay
  })
  observeEvent(input$typeG,{
    x$Gestion$typeG <- match(input$typeG, c("+", "x")) -1
  })
  observeEvent(input$upd,{
    x$Gestion$upd <- match(input$upd, c("Yes", "No"))
  })
  observeEvent(input$sup,{
    x$Gestion$sup <- input$sup
  })
  observeEvent(input$inf,{
    x$Gestion$inf <- input$inf
  })
  # tables
  output$tac <- renderRHandsontable({
    DF = tac()
    if (!is.null(DF)){
      x <- rhandsontable(DF, stretchH = "all")
      x <- hot_validate_numeric(
        x, cols = 1:ncol(DF), min = -1000, max = 1000, allowInvalid = TRUE
      )
    }
    x
  })
  observeEvent(input$tac, {
    x$Gestion$tac <- Gestion_tab$TAC[1,]
    x$Gestion$fbar <- Gestion_tab$TAC[2,]
  })

  output$eff <- renderRHandsontable({
    DF = maxeff()
    if (!is.null(DF)){
      x <- rhandsontable(DF, stretchH = "all", rowHeaderWidth = 150)
      x <- hot_validate_numeric(
        x, cols = 1:ncol(DF), min = -1000, max = 1000, allowInvalid = TRUE
      )
    }
    x
  })
  observeEvent(input$eff, {
    x$Gestion$effSup <- Gestion_tab$eff
  })

  output$mfm <- renderRHandsontable({
    DF = mfm()
    if (!is.null(DF)){
      x <- rhandsontable(DF, stretchH = "all", rowHeaderWidth = 150)
      x <- hot_validate_numeric(
        x, cols = 1:ncol(DF), min = -1000, max = 1000, allowInvalid = TRUE
      )
    }
    x
  })
  observeEvent(input$mfm, {
    x$Gestion$mfm <- Gestion_tab$mfm
  })
  # EoF Gestion

  ## Economic ####
  observeEvent(input$perscCalc, {
    x$Eco$perscCalc <- as.numeric(input$perscCalc)
  })
  observeEvent(input$dr, {
    x$Eco$dr <- input$dr
  })
  # EoF Economic

  ## Scenario ####
  observeEvent(input$scenar, {
    toggleState(id = "scen_var", condition = input$scenar)
    x$Scenario$active <- as.numeric(input$scenar)
  })
  observeEvent(input$scen_var, {
    x$Scenario$SELECTscen <- match(input$scen_var, x$Scenario$ALLscenario)
  })
  # EoF Scenario

  observe({
    # TODO : maybe only modify input at the very end.
    x$input@arguments$Recruitment <- x$Recruitment
    x$input@arguments$Replicates <- x$Replicates
    x$input@arguments$Scenario <- x$Scenario
    x$input@arguments$Gestion <- x$Gestion
    x$input@arguments$Eco <- x$Eco
    # print(summary(x$input))
  })

  # End of app ####
  observeEvent(input$done, {
    # cat("You're really not going to like it\n")
    # returnValue <- 42
    returnValue <- x$input
    stopApp(returnValue)
  })
  observeEvent(input$cancel, {
    # cat("Life? Don't talk to me about life.\n")
    stopApp()
  })
}


#' @noRd
#' @param object \code{\link[IAM]{iamInput-class}} or
#' \code{\link[IAM]{iamArgs-class}} object.
#' @param AllVarRep list of variable Outputrep can produce.
IAM_arg_app <- function(object, AllVarRep) {

  res <- with_golem_options(
    app = shinyApp(ui = args_app_ui, server = args_app_server),
    golem_opts = list(input = object, AllVarRep = AllVarRep,
                      spe = object@specific, arg = object@arguments,
                      tac = rbind(tac = object@arguments$Gestion$tac,
                                  fbar = object@arguments$Gestion$fbar),
                      maxeff = object@arguments$Gestion$effSup,
                      mfm =  object@arguments$Gestion$mfm)
  )
  if(is.null(res)){
    return(invisible())
  } else {
    return(res)
  }
}




#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#M?thodes de constitution des objets 'iamArgs' par interface graphique

#Arguments...

#' 'iamArgs' objects creator by GUI
#'
#' @param object \code{\link[IAM]{iamInput-class}} or
#' \code{\link[IAM]{iamArgs-class}} object.
#' @param specific this is for the specific file, a deprecated use of
#' this function.
#' @param ... desc parameter described below.
#'
#'  @return \code{\link[IAM]{iamArgs-class}} object after modification by a GUI.
#'
#' @docType methods
#' @name IAM.args-methods
#' @aliases IAM.args
#'
#' @examples
#' \dontrun{
#' data(IAM_input_2009)
#' # example not run because depend on shiny app.
#' # Initiate argument from input.
#' IAM_argum_2009 <- IAM.args(IAM_input_2009)
#' summary(IAM_argum_2009)
#' # Modify argument from argument.
#' IAM_argum_2009 <- IAM.args(IAM_argum_2009)
#' summary(IAM_argum_2009)
#' }
#'
#' @export
setGeneric("IAM.args", function(object, specific,  ...){ # Generic method ####
	standardGeneric("IAM.args")
	}
)


#' Etape d'initialisation
#' @param desc Object descriptor (default value : \code{as.character(NA)}).
#' If not provided, copied the description slot of object.
#' @rdname IAM.args-methods
setMethod("IAM.args", signature("iamInput","missing"),
          function(object, desc=as.character(NA), ...){ # input method ####

  if(is.null(desc)){ desc <- object@desc }

  args <- IAM.input2args(object, desc = desc)
  IAM.args(object = args, desc = desc)

})

#' Etape de modification
#' @rdname IAM.args-methods
setMethod("IAM.args", signature("iamArgs","missing"),
          function(object, desc=as.character(NA), ...){ # argum method ####
            if(is.null(desc)){ desc <- object@desc }

            AllVarRep = c(
              "B", "SSB", "Ctot", "Ytot", "Yfmi", "Ffmi", "Zeit", "Fbar", "Foth",
              "mu_nbds", "mu_nbv", "N", "Ystat", "Lstat", "Dstat", "Eff",
              "GVL_fme", "StatGVL_fme", "GVLtot_fm", "GVLav_f", "vcst_fm", "vcst_f",
              "rtbs_f", "gp_f", "ps_f", "gcf_f", "gva_f", "cs_f", "sts_f", "rtbsAct_f",
              "csAct_f", "gvaAct_f", "gcfAct_f", "psAct_f", "stsAct_f", "ccwCr_f",
              "GVLtot_f", "wagen_f", "L_efmit", "D_efmit",
              "Fr_fmi", "C_efmit", "P", "Pstat"
            )

            tac_dimcst <- attributes(object@arguments$Gestion$tac)$DimCst
            fbar_dimcst <- attributes(object@arguments$Gestion$fbar)$DimCst
            effsup_dimcst <- attributes(object@arguments$Gestion$effSup)$DimCst
            mfm_dimcst <- attributes(object@arguments$Gestion$mfm)$DimCst

            res <- IAM_arg_app(object = object, AllVarRep = AllVarRep)

            attributes(res@arguments$Gestion$tac)$DimCst <- tac_dimcst
            attributes(res@arguments$Gestion$fbar)$DimCst <- fbar_dimcst
            attributes(res@arguments$Gestion$effSup)$DimCst <- effsup_dimcst
            attributes(res@arguments$Gestion$mfm)$DimCst <- mfm_dimcst

            return(res)
          })


#' Etape d'initialisation from txt file.
#' @rdname IAM.args-methods
#' @importFrom methods new
setMethod("IAM.args", signature("character","character"),
          function(object, specific, desc=as.character(NA), ...){ # depr method ####

  if (substring(object,nchar(object)-3,nchar(object))!=".txt") stop("'object' must be a .txt file!!")
  if (substring(specific,nchar(specific)-3,nchar(specific))!=".txt") stop("'specific' must be a .txt file!!")

  stop("This depend on deprecated C++ function 'Fun' ")
  # specif <- .Call("Fun",normalizePath(specific),NULL)
  # argum <- .Call("Fun",normalizePath(object),specif)

  # return(new("iamArgs", desc=desc, arguments=argum, specific=specif))

})

#:::::::::::::
#Examples
#:::::::::::::
##importation des arguments
#impArg <- IAM.args("C:/Documents and Settings/mmerzere/Bureau/COST_R/IAMwd/expArgs/args.txt",
#                   "C:/Documents and Settings/mmerzere/Bureau/COST_R/IAMwd/expArgs/specific.txt",desc="My args")
#
