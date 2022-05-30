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
#' data(IAM_argum_2009)
#' IAM_argum_2009@arguments$Eco$dr
#' IAM_argum_2009@arguments$Eco$perscCalc
#' IAM_argum_2009 <- IAM.editArgs_Eco(IAM_argum_2009, dr = 0.06, perscCalc = 0)
#' IAM_argum_2009@arguments$Eco$dr
#' IAM_argum_2009@arguments$Eco$perscCalc
#'
#' @name IAM.editArgs_Eco
#' @export
setGeneric("IAM.editArgs_Eco",
           function(object, ...){ ## Generic editArgs ####
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
  function(object, dr = NULL, perscCalc = NULL) { ## iamArgs editArgs ####

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
#' data(IAM_argum_2009)
#' IAM.args_scenar(IAM_argum_2009)
#'
#'
#' @name IAM.args_scenar
#' @export
setGeneric("IAM.args_scenar",
           function(object){ ## generic arg scenar ####
  standardGeneric("IAM.args_scenar")
}
)

#' @aliases IAM.args_scenar,iamArgs-method
#' @rdname IAM.args_scenar
#' @export
setMethod(
  f = "IAM.args_scenar", signature("iamArgs"),
  function(object) { ## iamArgs arg scenar ####
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
#' data(IAM_argum_2009)
#' IAM.args_scenar(IAM_argum_2009)
#'
#' # Activate scenario "Sc_saisine"
#' IAM_argum_2009 <- IAM.editArgs_Scenar(IAM_argum_2009, selected="Sc_saisine")
#' IAM.args_scenar(IAM_argum_2009)
#'
#' # Desactivate scenario
#' IAM_argum_2009 <- IAM.editArgs_Scenar(IAM_argum_2009, selected = NULL)
#' IAM.args_scenar(IAM_argum_2009)
#'
#' # Activate scenario "Sc_saisine"
#' IAM_argum_2009 <- IAM.editArgs_Scenar(IAM_argum_2009, selected = 1)
#' IAM.args_scenar(IAM_argum_2009)
#'
#'
#' @name IAM.editArgs_Scenar
#' @export
setGeneric("IAM.editArgs_Scenar",
           function(object, ...){ ## generic editArgs ####
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
  function(object, selected = NULL) { ## iamArgs editArgs ####
    scenar <- IAM.args_scenar(object)
    scenarii <- scenar$ALLscenario

    if(!is.null(selected)){
      if(is.numeric(selected) && selected %in% seq_along(scenarii)){
        scenar$SELECTscen <- selected
        scenar$active <- 1L
      } else if(is.character(selected) && selected %in% scenarii){
        scenar$SELECTscen <- match(selected, scenarii)
        scenar$active <- 1L
      } else {
        stop(paste("selected must be either a numeric between 1 and the",
                   " number of scenarii or the name of a scenario"))
      }
    } else {
      scenar$active <- 0L
    }

    object@arguments$Scenario <- scenar
    return(object)
  }
)

# Gestion ####
#' Modify Gestion arguments
#'
#' @param object a \code{\link[IAM]{iamArgs-class}} object
#' @param ... selected argument
#'
#' @include 0_input_objects.r
#'
#' @author Maxime Jaunatre
#'
#' @examples
#' data(IAM_argum_2009)
#' IAM.args_scenar(IAM_argum_2009)
#'
#' # Activate scenario "Sc_saisine"
#' IAM_argum_2009 <- IAM.editArgs_Gest(IAM_argum_2009, active= TRUE)
#' summary(IAM_argum_2009)
#'
#'
#' @name IAM.editArgs_Gest
#' @export
setGeneric("IAM.editArgs_Gest",
           function(object, ...){ ## generic editArgs ####
             standardGeneric("IAM.editArgs_Gest")
           }
)

#' @param active Activation of the module. TRUE or FALSE input,
#' but accept 1 or 0 values (any numeric will be interpreted as logical)
#' @param control option to choose if the management will control effort
#' reduction by reducing the number of vessels or number of trips per vessels.
#' @param target management aim to reduce effort to reach TAC or Fbar,
#' or TAC then Fbar.
#' @param espece species to target for management.
#' @param delay option to delay the management for n years. If the delay is
#' larger than the year in model, the management will have no effect.
#' @param type option to choose if effort reduction is additive ("+") or
#' multiplicative ("x").
#' @param update option to choose is effort reduction is from t0 or t-1
#' iteration
#' @param bounds bounds for target reach at each model iteration.
#' num of length 2.
#' @param tac TAC values for each model iteration. Can be NA.
#' @param fbar Fbar values for each model iteration. Can be NA.
#' @param effSup Maximum effort per fleet and year. To limit nbds.
#' @param mfm Ponderation fleet/metier for effort reduction.
#'
#'
#' @rdname IAM.editArgs_Gest
#' @export
setMethod(
  f = "IAM.editArgs_Gest", signature("iamArgs"),
  function(object, active = NULL,
           control = NULL,
           target = NULL, espece = NULL,
           delay = NULL, type = NULL, update = NULL, bounds = NULL,
           tac = NULL, fbar = NULL, effSup = NULL, mfm = NULL) { ## iamArgs editArgs ####

    gest <- object@arguments$Gestion

    if(!is.null(active)){
      if(is.numeric(active)){
        active <- as.logical(active)
      }
      if(!is.logical(active) | length(active) > 1 ) {
        stop(paste("active must be logical (or numeric, FALSE being 0 and TRUE",
                   "being 1) single value"), call. = interactive())
      }
      gest$active <- as.numeric(active)
    }

    if(!is.null(control)){
      tryCatch(
        control <- match.arg(control, c("Nb vessels", "Nb trips")),
        error = function(c) stop(paste("control must be either 'Nb vessels' or",
                                       "'Nb trips'", call. = FALSE))
      )
      gest$control <- control
    }

    if(!is.null(target)){
      tryCatch(
        target <- match.arg(target, c("TAC", "Fbar", "TAC->Fbar")),
        error = function(c) stop(paste("target must be either 'TAC', 'Fbar' or",
                                       "'TAC->Fbar'", call. = FALSE))
      )
      gest$target <- target
    }

    if(!is.null(espece)){
      espece <- match.arg(espece, object@specific$Species) # TODO :All species or only XSA ?
      gest$espece <- espece
    }

    if(!is.null(delay)){
      if(!is.numeric(delay) | length(delay) > 1 | is.na(delay)| delay < 0 ) {
        stop("delay must be a positive numeric single value",
             call. = interactive())
      }
      if(delay > object@specific$NbSteps){
        warning(paste("delay value is greater than number of steps and thus ",
                      "scenario will never be triggered"), call. = interactive())
      }
      gest$delay <- as.integer(delay)
    }

    if(!is.null(type)){
      tryCatch(
        type <- match.arg(type, c("+", "x")),
        error = function(c) stop("type must be either '+' or 'x'", call. = FALSE)
      )
      type <- match(type, c("+", "x")) - 1
      gest$typeG <- type
    }

    if(!is.null(update)){
      warning("update argument for Gestion module is deprecated !")
      if(is.numeric(update)){
        update <- as.logical(update)
      }
      if(!is.logical(update) | length(update) > 1 ) {
        stop(paste("update must be logical (or numeric, FALSE being 0 and TRUE",
                   "being 1) single value"), call. = interactive())
      }
      gest$upd <- -as.numeric(update) + 2
    }

    if(!is.null(bounds)){
      if(!is.numeric(bounds) | length(bounds) != 2 ) {
        stop(paste("bounds must be 2 numeric values"), call. = interactive())
      }
      # bounds will be sorted
      bounds <- sort(bounds)

      gest$sup <- bounds[2]
      gest$inf <- bounds[1]
    }

    if(!is.null(tac)){
      n <- length(gest$tac)
      if(any(!is.numeric(tac)) | length(tac) != n |
         any(tac[!is.na(tac)] < 0) ) {
        stop(sprintf("tac must be %d positive numeric values", n),
             call. = interactive())
      }
      gest$tac[] <- tac
    }

    if(!is.null(fbar)){
      n <- length(gest$fbar)
      if(any(!is.numeric(fbar)) | length(fbar) != n |
         any(fbar[!is.na(fbar)] < 0) ) {
        stop(sprintf("fbar must be %d positive numeric values", n),
             call. = interactive())
      }
      gest$fbar[] <- fbar
    }

    if(!is.null(effSup)){
      n <- dim(gest$effSup)
      if(any(!is.numeric(effSup)) | any(dim(effSup) != n) |
         any(effSup[!is.na(effSup)] < 0) ) {
        stop(sprintf("effSup must be a %d by %d positive numeric matrix",
                     n[1], n[2]),
             call. = interactive())
      }
      gest$effSup[] <- effSup
    }

    if(!is.null(mfm)){
      n <- dim(gest$mfm)
      if(any(!is.numeric(mfm)) | any(dim(mfm) != n)  ) {
        stop(sprintf("mfm must be a %d by %d numeric matrix",
                     n[1], n[2]),
             call. = interactive())
      }
      gest$mfm[] <- mfm
    }


    object@arguments$Gestion <- gest
    return(object)
  }
)


# Replicates ####
#' Modify Replicates arguments
#'
#' @param object a \code{\link[IAM]{iamArgs-class}} object
#' @param ... selected argument
#'
#' @include 0_input_objects.r
#'
#' @author Maxime Jaunatre
#'
#' @examples
#' data(IAM_argum_2009)
#'
#' # Replicate 100 times
#' IAM_argum_2009 <- IAM.editArgs_Rep(IAM_argum_2009, n = 100)
#' IAM.args_scenar(IAM_argum_2009)
#'
#' # Desactivate Replicates
#' IAM_argum_2009 <- IAM.editArgs_Rep(IAM_argum_2009, n = NULL)
#' IAM.args_scenar(IAM_argum_2009)
#'
#' @name IAM.editArgs_Rep
#' @export
setGeneric("IAM.editArgs_Rep",
           function(object, ...){ ## generic editArgs ####
             standardGeneric("IAM.editArgs_Rep")
           }
)

#' @param n Number of replicates
#' NULL will desactivate the Replicates module
#'
#' @rdname IAM.editArgs_Rep
#' @export
setMethod(
  f = "IAM.editArgs_Rep", signature("iamArgs"),
  function(object, n = NULL) { ## iamArgs editArgs ####
    rep <- object@arguments$Replicates

    if(!is.null(n)){
      if(is.numeric(n) && length(n) == 1 && n > 0){
        rep$nbIter <- n
        rep$active <- 1L
      } else {
        stop("n must be a single positive numeric value")
      }
    } else {
      rep$active <- 0L
    }

    object@arguments$Replicates <- rep
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
#' data("IAM_argum_2009")
#' IAM_argum_2009@arguments$Gestion$active <- 1
#' IAM_argum_2009@arguments$Replicates$active <- 1
#' IAM_argum_2009@arguments$Scenario$ALLscenario <- c("Scenario1", "Scenario2")
#' IAM_argum_2009 <- IAM.editArgs_Scenar(IAM_argum_2009, selected = 1)
#' summary(IAM_argum_2009)
#'
#' @export
setMethod(
  f = "summary", signature("iamArgs"),
  function(object, ...) { ## iamArgs summary ####
    spe <- object@specific
    arg <- object@arguments

    if (any(spe$Q == 1)) {
      listType <- list(
        modSRactive=0,typeMODsr="Mean",parAmodSR=0,parBmodSR=0,parCmodSR=0,
        wnNOISEmodSR=0,noiseTypeSR=1,simuSTOCHactive=0,typeSIMUstoch=1
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
        rep("=",87),
        "\n  SR module  |               Stock Recruitment             |",
        "      Noise       | Proba |\n",
        rep("-",87),
        "\n    Species  |    function  :  param A ; param B ; param C |",
        "  Type :    sd    |  Type |\n",
        sprintf(" %5s (%3s) | %13s %.3e  %.2e  %.2e |  %4s | %.2e |   %1s   |\n",
                spe$Species,
                ifelse(spe$Q == 0, "XSA", "SS3"),
                ifelse(SRactif == 1, unlist(lapply(arg$Recruitment, "[[", "typeMODsr")),
                       "not activated" ),
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
        rep("-",87),"\n",
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
#' data("IAM_input_2009")
#' summary(IAM_input_2009)
#'
#' @export
setMethod(
  f = "summary", signature("iamInput"),
  function(object, ...) { ## iamInput summary ####
    spe <- object@specific

    cat(object@desc, "(IAM input) :\n")

    cat(sprintf(
      "Simulation of %d dynamic species, %d static species and %d fleet\n",
      length(spe$Species), length(spe$StaticSpp), length(spe$Fleet)
    ))
    cat(sprintf("Simulation start in %d and end in %d (%d steps)\n",
                spe$t_init, tail(spe$times,1), spe$NbSteps ))

    cat("\n",rep("-",41), "\nDynamic Species | Model |     Ages (nb) |\n", sep = "")
    cat(sprintf("%15s | %5s | %1s to %3s (%2d) |\n", spe$Species,
                ifelse(spe$Q == 0, "XSA", "SS3"),
                unlist(lapply(spe$Ages, function(x) head(x,1))),
                unlist(lapply(spe$Ages, function(x) tail(x, 1))),
                unlist(lapply(spe$Ages, function(x) length(x)))), sep = "")

    cat(rep("-",41),"\n", sep = "")

    if(length(spe$Fleet) < 15){
      cat("                  Fleet |    nbv   |\n", sprintf(
        "%22s |  %4d    |\n", spe$Fleet, object@input$Fleet$nbv_f
      ))
    } else {
      cat(" Too mant fleet to summarise ( > 15 )\n")
    }


  }
)





