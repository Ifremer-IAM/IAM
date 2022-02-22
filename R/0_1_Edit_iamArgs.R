# Economic ####
#' Modify Economic arguments
#'
#' @param object a \code{\link[IAM]{iamArgs-class}} object
#' @param ... dr and persCalc arguments
#'
#' @details
#' The GUI provide many other arguments for economic parameters but they are
#' deprecated and it is useless to edit them
#'
#' @include 0_input_objects.r
#'
#' @author Maxime Jaunatre
#'
#' @examples
#' data(IAM_argum_1984)
#' IAM_argum_1984@arguments$Eco$dr
#' IAM_argum_1984@arguments$Eco$perscCalc
#' IAM_argum_1984 <- IAM.editArgs_Eco(IAM_argum_1984, dr = 0.06, perscCalc = 0)
#' IAM_argum_1984@arguments$Eco$dr
#' IAM_argum_1984@arguments$Eco$perscCalc
#'
#' @name IAM.editArgs_Eco
#' @export
setGeneric("IAM.editArgs_Eco", function(object, ...){
  standardGeneric("IAM.editArgs_Eco")
}
)

#' @param dr actualisation rate : numeric value in decimal
#' (percentage. For 4\%, dr must be 0.04)
#' @param perscCalc personnel cost : computation mode
#' \describe{
#' \item{0}{- constant variable}
#' \item{1}{- linear relation with rtbs; cshr available}
#' \item{2}{- linear relation with rtbs; cshr non available}
#' }
#'
#' @rdname IAM.editArgs_Eco
#' @export
setMethod(
  f = "IAM.editArgs_Eco", signature("iamArgs"),
  function(object, dr = NULL, perscCalc = NULL) {

    eco <- object@arguments$Eco

    if(!is.null(dr)){
      if(!is.numeric(dr) | length(dr) > 1 ) {
        stop("dr must be a numeric single value", call. = interactive())
      }
      eco$dr <- dr
    }

    if(!is.null(perscCalc)){
      if(!is.numeric(perscCalc) | length(perscCalc) > 1 ) {
        stop("persCalc must be a numeric single value", call. = interactive())
      }
      eco$perscCalc <- as.integer(perscCalc)
    }
    object@arguments$Eco <- eco
    return(object)
  }
)



# Scenario ####
#' Get Scenario arguments
#'
#' Return the scenario arguments.
#'
#' @param object a \code{\link[IAM]{iamArgs-class}} object
#'
#' @include 0_input_objects.r
#'
#' @author Maxime Jaunatre
#'
#' @examples
#' data(IAM_argum_1984)
#' IAM.args_scenar(IAM_argum_1984)
#'
#'
#' @name IAM.args_scenar
#' @export
setGeneric("IAM.args_scenar", function(object){
  standardGeneric("IAM.args_scenar")
}
)

#' @aliases IAM.args_scenar,iamArgs-method
#' @rdname IAM.args_scenar
#' @export
setMethod(
  f = "IAM.args_scenar", signature("iamArgs"),
  function(object) {
    object@arguments$Scenario
  }
)



#' Modify Scenario arguments
#'
#' @param object a \code{\link[IAM]{iamArgs-class}} object
#' @param ... selected argument
#'
#' @include 0_input_objects.r
#'
#' @author Maxime Jaunatre
#'
#' @examples
#' data(IAM_argum_1984)
#' IAM.args_scenar(IAM_argum_1984)
#'
#' # Activate scenario "Sc_saisine"
#' IAM_argum_1984 <- IAM.editArgs_Scenar(IAM_argum_1984, selected="Sc_saisine")
#' IAM.args_scenar(IAM_argum_1984)
#'
#' # Desactivate scenario
#' IAM_argum_1984 <- IAM.editArgs_Scenar(IAM_argum_1984, selected = NULL)
#' IAM.args_scenar(IAM_argum_1984)
#'
#' # Activate scenario "Sc_saisine"
#' IAM_argum_1984 <- IAM.editArgs_Scenar(IAM_argum_1984, selected = 1)
#' IAM.args_scenar(IAM_argum_1984)
#'
#'
#' @name IAM.editArgs_Scenar
#' @export
setGeneric("IAM.editArgs_Scenar", function(object, ...){
  standardGeneric("IAM.editArgs_Scenar")
}
)

#' @param selected name or number of scenario selected.
#' NULL will desactivate the scenario module
#'
#' @rdname IAM.editArgs_Scenar
#' @export
setMethod(
  f = "IAM.editArgs_Scenar", signature("iamArgs"),
  function(object, selected = NULL) {
    scenar <- IAM.args_scenar(object)
    scenarii <- scenar$ALLscenario

    if(!is.null(selected)){
      if(is.numeric(selected) && selected %in% seq_along(scenarii)){
        scenar$SELECTscen <- selected
        scenar$active <- 1
      } else if(is.character(selected) && selected %in% scenarii){
        scenar$SELECTscen <- match(selected, scenarii)
        scenar$active <- 1
      } else {
        stop(paste("selected must be either a numeric between 1 and the",
                   " number of scenarii or the name of a scenario"))
      }
    } else {
      scenar$active <- 0
    }

    object@arguments$Scenario <- scenar
    return(object)
  }
)


# Summary ####
#' summary method for iamInput class
#'
#' @param object a \code{\link[IAM]{iamArgs-class}} object
#' @param ... Unused argument
#'
#' @details
#' The GUI provide many other arguments for economic parameters but they are
#' deprecated and it is useless to show them
#'
#' @include 0_input_objects.r
#' @importFrom utils tail
#'
#' @author Maxime Jaunatre
#'
#' @examples
#' data("IAM_argum_1984")
#' IAM_argum_1984@arguments$Gestion$active <- 1
#' IAM_argum_1984@arguments$Replicates$active <- 1
#' IAM_argum_1984@arguments$Scenario$ALLscenario <- c("Scenario1", "Scenario2")
#' IAM_argum_1984 <- IAM.editArgs_Scenar(IAM_argum_1984, selected = 1)
#' summary(IAM_argum_1984)
#'
#' @export
setMethod(
  f = "summary", signature("iamArgs"),
  function(object, ...) {
    spe <- object@specific
    arg <- object@arguments

    if (any(spe$Q == 1)) {
      listType <- list(
        modSRactive=0,typeMODsr="Mean",parAmodSR=NA,parBmodSR=NA,parCmodSR=NA,
        wnNOISEmodSR=NA,noiseTypeSR=1,simuSTOCHactive=0,typeSIMUstoch=1
        )
      namQ <- names(spe$Q)[spe$Q == 1]
      lll <- lapply(namQ,function(z) return(listType))
      names(lll) <- namQ
      newL <- c(arg$Recruitment,lll)
      arg$Recruitment <- newL[spe$Species]
    } # Eof SS3 SR

    SRactif <- unlist(lapply(arg$Recruitment, "[[", "modSRactive"))
    Stockactif <- unlist(lapply(arg$Recruitment, "[[", "simuSTOCHactive"))
    TypNoise <- unlist(lapply(arg$Recruitment, "[[", "noiseTypeSR"))

    cat(object@desc, "(IAM argument) :\n")
    cat(sprintf(
      "Simulation of %d dynamic species, %d static species and %d fleet\n",
      length(spe$Species), length(spe$StaticSpp), length(spe$Fleet)
    ))
    cat(sprintf("Simulation start in %d and end in %d (%d steps)\n",
                spe$t_init, tail(spe$times,1), spe$NbSteps ))

    cat("\n",
        rep("=",83),
        "\n  SR module  |               Stock Recruitment             |",
        "    Noise     | Proba |\n",
        rep("-",83),
        "\n    Species  |    function  :  param A ; param B ; param C |",
        "  Type :  sd  |  Type |\n",
        sprintf(" %5s (%3s) | %13s %9d  %8d  %8d |  %4s |%.2g |   %1s   |\n",
                spe$Species,
                ifelse(spe$Q == 0, "XSA", "SS3"),
                ifelse(SRactif == 1, unlist(lapply(arg$Recruitment, "[[", "typeMODsr")),
                       "............." ),
                unlist(lapply(arg$Recruitment, "[[", "parAmodSR")), # A
                unlist(lapply(arg$Recruitment, "[[", "parBmodSR")), # B
                unlist(lapply(arg$Recruitment, "[[", "parCmodSR")), # C
                ifelse(TypNoise == 1, "Norm", "LogN" ), # Type Noise
                unlist(lapply(arg$Recruitment, "[[", "wnNOISEmodSR")), # Noise value
                ifelse(
                  Stockactif == 1,
                  as.character(unlist(lapply(arg$Recruitment, "[[", "typeSIMUstoch"))), "."
                  )
        ),
        rep("-",83),"\n",
        sep = "")

    if(arg$Gestion$active == 1){
      cat("\n",
          rep("=",77),
          "\n                     Gestion Module                      ",
          "                   | \n",
          rep("-",77),
          "\n    Species  |   control  |    target | delay | typeG | update |",
          "   bounds   |\n",
          sprintf(
            " %5s (%3s) | %10s | %9s |  %3d  |   %1s   |   %3s  | %3d : %4d |\n",
            arg$Gestion$espece,
            ifelse(spe$Q[arg$Gestion$espece] == 0, "XSA", "SS3"),
            arg$Gestion$control,
            arg$Gestion$target,
            arg$Gestion$delay,
            ifelse(arg$Gestion$typeG == 0, "+", "x" ),
            ifelse(arg$Gestion$upd == 1, "yes", "no" ),
            arg$Gestion$sup, arg$Gestion$inf
          ),
          rep("-",77),"\n",
          sep = "")
      cat(
        "   TAC :", arg$Gestion$tac, "\n",
        " FBAR :", arg$Gestion$fbar,"\n",
        sep = " ")
      cat(rep("-",76),"\n", sep = "")
    } else {
      cat("\n The Gestion module is not active.\n")
    }

    # Eco and replicates
    cat("\n",
        rep("=", 60), "\n",
        sprintf(
          "Economic : PerscCalc = %1d ; dr = %0.3f ",
          arg$Eco$perscCalc, arg$Eco$dr
        ),
        sep = ""
    )
    if (arg$Replicates$active == 1) {
      cat(
        sprintf("| Replicates : %5d |\n", arg$Replicates$nbIter),
        sep = ""
      )
    } else {
      cat("| No replicates      |\n")
    }
    cat(rep("-", 60), "\n", sep = "")

    # Ajouter les scenario
    if(arg$Scenario$active == 1){
      cat("\n------------------------")
      cat("\n Scenario (selected:*) |")
      for(sce in seq_along(arg$Scenario$ALLscenario)){
        cat(sprintf(
          "\n %19s %1s |",
          arg$Scenario$ALLscenario[sce],
          ifelse(sce == arg$Scenario$SELECTscen, "*", "")), sep = "")
      }
      cat("\n------------------------\n")
    } else {
      cat("\n The Scenario module is not active.\n")
    }

  }
)



#' summary method for iamInput class
#'
#' @param object a \code{\link[IAM]{iamInput-class}} object
#' @param ... Unused argument
#'
#' @include 0_input_objects.r
#'
#' @author Maxime Jaunatre
#'
#' @importFrom utils tail head
#'
#' @examples
#' data("IAM_input_1984")
#' summary(IAM_input_1984)
#'
#' @export
setMethod(
  f = "summary", signature("iamInput"),
  function(object, ...) {
    spe <- object@specific

    cat(object@desc, "(IAM input) :\n")

    cat(sprintf(
      "Simulation of %d dynamic species, %d static species and %d fleet\n",
      length(spe$Species), length(spe$StaticSpp), length(spe$Fleet)
    ))
    cat(sprintf("Simulation start in %d and end in %d (%d steps)\n",
                spe$t_init, tail(spe$times,1), spe$NbSteps ))

    cat("\n",rep("-",36), "\nDynamic Species | Model |     Ages |\n", sep = "")
    cat(sprintf("%15s | %5s | %1s to %3s |\n", spe$Species,
                ifelse(spe$Q == 0, "XSA", "SS3"),
                unlist(lapply(spe$Ages, function(x) head(x,1))),
                unlist(lapply(spe$Ages, function(x) tail(x, 1))) ), sep = "")

    cat(rep("-",36),"\n", sep = "")
    cat("                  Fleet |    nbv   |\n", sprintf(
      "%22s |  %4d    |\n", spe$Fleet, object@input$Fleet$nbv_f
    ))

  }
)


