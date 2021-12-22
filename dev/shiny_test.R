library(shiny)
library(shiny.i18n)
library(shinyWidgets)
library(shinythemes)
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
            ),
            div( # Scenar ####
                 class = "option-group",
                 checkboxInput("scenar", "Scenario", value = FALSE),
                 selectInput("scen_var", "Output variables",
                             get_golem_options("input")@arguments$Scenario$ALLscenario,
                             selectize = FALSE
                 )
            )
          ),
          column(
            width = 6,
            div( # Economic panel ####
                 class = "option-group",
                 h3("Economic"),
                 radioGroupButtons(
                   inputId = "Eco_typ", choices = c("Complete","DCF"),
                   justified = TRUE
                 ),
                 awesomeRadio(
                   inputId = "pers", label = "perscCalc", choices = c("0", "1", "2", "3", "4"),
                   selected = "0", inline = TRUE
                 ),
                 numericInput("disc_rate", label = "Discound rate", value = 0, step = 0.01)
                 # )
            ),
            div( # Iterative panel ####
                 class = "option-group",
                 checkboxInput("iter", "Iterative", value = FALSE),
                 numericInput("niter", label = "Number of iteration", value = 500),
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
        actionButton("done", "Done"),
        actionButton("cancel", "Cancel")
      )
    )
  )
}

# Define server function
app_server <- function(input, output, session) {
  sp <- get_golem_options("input")@specific$Species
  spDyn <- sp[get_golem_options("input")@specific$Q == 0]
  for( i in spDyn){
    appendTab(inputId = "tabs",
              tabPanel(i, mod_spInput(i)) # TODO : need to reload interface with default once tab is changed
    )
  }
  updateTabsetPanel(session, "tabs", selected = unname(sp[1]))


  observeEvent(input$iter, {
    toggleState(id = "niter")
    toggleState(id = "var_rep")
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

  # End of app ####
  observeEvent(input$done, {
    cat("You're really not going to like it\n")
    returnValue <- 42
    stopApp(returnValue)
    # TODO : return args and modify all it's a values according to what's in shiny
  })
  observeEvent(input$cancel, {
    cat("Life? Don't talk to me about life.\n")
    stopApp()
  })
}

run_app <- function(object, AllVarRep) {
  res <- with_golem_options(
    app = shinyApp(ui = app_ui, server = app_server),
    golem_opts = list(input = object, AllVarRep = AllVarRep)
  )
  if(is.null(res)){
    return(invisible())
  } else {
    return(res)
  }
}

# Plotting Rec ####
# 1 -> recrutement constant moyen (rec~a)
# 2 -> Hockey stick (rec ~ (si (ssb<=b) a*ssb sinon a*b))
# 3 -> Beverton & Holt (rec ~ a*ssb/(b+ssb))
# 4 -> Ricker (rec ~ a*ssb*exp(-b*ssb))
# 5 -> Shepherd (rec ~ a*ssb/(1+ (ssb/b)^c))
# 6 -> Hockey Stick Quadratic (rec ~ (si (ssb<=b*(1-c)) a*ssb ; si (b*(1-c)<ssb<b*(1+c)) a*(ssb-((ssb-b*(1-c))^2)/(4*b*c)) ; sinon a*b))
# 7 -> Hockey Stick Smooth (rec ~ a*(ssb+sqrt(b^2+g)-sqrt((ssb-b)^2+g)), avec g=0.001 )


# TODO : what could be very nice is to initialise iam.Args
# with the gui default and then just write a function that modify it.
# So we could call it very quickly in script for automated test
setMethod("IAM.input2args", signature("iamInput","missing"),function(object, desc=as.character(NA), ...){

  if(is.null(desc)){ desc <- object@desc }
  # init the arg object with shiny default
  ## Create argum ####
  # TODO : this is where to modify the default for the GUI now !
  init_recru <- function(name, input){
    list(modSRactive = 1,
         typeMODsr = "Mean",
         parAmodSR = unname(input@input[[name]]$N_it0[1]),
         parBmodSR = 0, parCmodSR = 0, wnNOISEmodSR = 0, noiseTypeSR = 1,
         simuSTOCHactive = 0, typeSIMUstoch = 1
    )
  }

  sp <- as.list(input@specific$Species) ; names(sp) <- sp
  Recruitement <- lapply(sp, function(x) init_recru(name = x, input=input))
  rm(sp)

  Replicates <- list(active =0, nbIter =500,
                     SELECTvar = ALLVarRep[c(1:20, 22:38, 43:44)] )

  Scenario <- list(active = 0, ALLscenario = names(input@scenario), SELECTscen = 1)

  Eco <- list(active = 0, type = 1,
              adj = 1, ue_choice = 1, oths = 0, othsFM = 0, # useless
              perscCalc = 0,
              report = 0, # useless
              dr = 0)

  init_gest <- function(input){
    years <- input@specific$times
    fleets <- input@specific$Fleet
    spp <- c(na.omit(input@specific$Species))
    sppStat <- c(na.omit(input@specific$StaticSpp))
    Espece <- c(na.omit(c(spp,sppStat)))[1]

    tacfbar <- matrix(as.numeric(NA),nrow=2,ncol=length(years),
                      dimnames=list(c("TAC","Fbar"),years))
    mfm <- input@input$Fleet$nbv_f_m ; mfm[!is.na(mfm)] <- 1

    TACOTH <- matrix(as.numeric(NA),nrow=length(sppStat)+length(spp),
                     ncol=length(years),dimnames=list(c(spp,sppStat),years))
    if (nrow(TACOTH)>1) {
      TACOTH <- as.data.frame(t(TACOTH[rownames(TACOTH)!=Espece,]))
      tabO <- lapply(TACOTH, function(x) {
        attributes(x)$DimCst <- as.integer(c(0,0,0,length(x))) ; return(x)})
    } else {
      tabO <- NULL
    }

    Gestion <- list(active = 0, control = "Nb vessels", target = "TAC",
                    espece = Espece, delay = 2, typeG = 0, upd = 1,
                    sup = 0, inf =0,
                    tac = tacfbar["TAC",],
                    fbar = tacfbar["Fbar",],
                    othSpSup = tabO,
                    effSup = matrix(as.numeric(NA),nrow=length(fleets),
                                    ncol=length(years),dimnames=list(fleets,years)),
                    mfm = mfm,
                    TACbyF = matrix(as.numeric(NA),nrow=length(fleets),
                                    ncol=length(years),dimnames=list(fleets,years)))

    attributes(Gestion$tac)$DimCst <- as.integer(c(0,0,0,length(Gestion$tac)))
    attributes(Gestion$fbar)$DimCst <- as.integer(c(0,0,0,length(Gestion$fbar)))
    attributes(Gestion$effSup)$DimCst <- as.integer(c(nrow(Gestion$effSup),0,0,
                                                      ncol(Gestion$effSup)))
    attributes(Gestion$mfm)$DimCst <- as.integer(c(dim(Gestion$mfm),0,0))
    attributes(Gestion$TACbyF)$DimCst <- as.integer(c(nrow(Gestion$TACbyF),0,0,
                                                      ncol(Gestion$TACbyF)))
    return(Gestion)
  }
  Gestion <- init_gest(input)

  argum <- list(Recruitment = Recruitement, Replicates = Replicates,
                Scenario = Scenario, Gestion = Gestion, Eco = Eco)
  ## Copy specific ####
  specif <- object@specific
  args <- new("iamArgs", desc = desc, arguments = argum, specific = specif)

  return(args)

})

setMethod("IAM.args", signature("iamInput","missing"),function(object, desc=as.character(NA), ...){

  if(is.null(desc)){ desc <- object@desc }

  args <- IAM.input2arg(object, desc = desc)
  IAM.args(object = args, desc = desc)

})

setMethod("IAM.args", signature("iamArgs","missing"),function(object, desc=as.character(NA), ...){
  if(is.null(desc)){ desc <- object@desc }

  ALLVarRep = c(
    "B", "SSB", "Ctot", "Ytot", "Yfmi", "Ffmi", "Zeit", "Fbar", "Foth", "mu_nbds", "mu_nbv", "N", "Ystat", "Lstat", "Dstat", "Eff",
    "GVL_fme", "StatGVL_fme", "GVLtot_fm", "GVLav_f", "vcst_fm", "vcst_f", "rtbs_f", "gp_f", "ps_f", "gcf_f", "gva_f", "cs_f", "sts_f", "rtbsAct_f",
    "csAct_f", "gvaAct_f", "gcfAct_f", "psAct_f", "stsAct_f", "ccwCr_f", "GVLtot_f", "wagen_f", "L_efmit", "D_efmit",
    "Fr_fmi", "C_efmit", "P", "Pstat"
  )
  return(run_app(object = object, ALLVarRep = ALLVarRep))
})


# just to test it for now
ALLVarRep = c(
  "B", "SSB", "Ctot", "Ytot", "Yfmi", "Ffmi", "Zeit", "Fbar", "Foth", "mu_nbds", "mu_nbv", "N", "Ystat", "Lstat", "Dstat", "Eff",
  "GVL_fme", "StatGVL_fme", "GVLtot_fm", "GVLav_f", "vcst_fm", "vcst_f", "rtbs_f", "gp_f", "ps_f", "gcf_f", "gva_f", "cs_f", "sts_f", "rtbsAct_f",
  "csAct_f", "gvaAct_f", "gcfAct_f", "psAct_f", "stsAct_f", "ccwCr_f", "GVLtot_f", "wagen_f", "L_efmit", "D_efmit",
  "Fr_fmi", "C_efmit", "P", "Pstat"
)

run_app(args, ALLVarRep)



