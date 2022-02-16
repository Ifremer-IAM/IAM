#' Modify Economic arguments
#'
#' @param object a \code{\link[IAM]{iamArgs-class}} object
#' @param dr actualisation rate : numeric value in decimal (percentage. For 4%, dr must be 0.04)
#' @param persCalc personnel cost : computation mode
#' \itemize{
#' \item{0 - constant variable}
#' \item{1 - linear relation with rtbs; cshr available}
#' \item{2 - linear relation with rtbs; cshr non available}
#' }
#'
#' @details
#' The GUI provide many other arguments for economic parameters but they are
#' deprecated and it is useless to edit them
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
#'
#' @name IAM.editArgs_Eco
#' @rdname IAM.editArgs_Eco
#' @export
setGeneric("IAM.editArgs_Eco", function(object, dr = NULL, perscCalc = NULL){
  standardGeneric("IAM.editArgs_Eco")
}
)

#' @name IAM.editArgs_Eco
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



#' Get Scenario arguments
#'
#' Return the scenario arguments.
#'
#' @param object a \code{\link[IAM]{iamArgs-class}} object
#'
#' @author Maxime Jaunatre
#'
#' @examples
#' data(IAM_argum_1984)
#' IAM.args_scenar(IAM_argum_1984)
#'
#'
#' @name IAM.args_scenar
#' @rdname IAM.args_scenar
#' @export
setGeneric("IAM.args_scenar", function(object){
  standardGeneric("IAM.args_scenar")
}
)

#' @name IAM.args_scenar
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
#' @param selected name or number of scenario selected.
#' NULL will desactivate the scenario module
#'
#' @details
#' The selection of
#'
#' @author Maxime Jaunatre
#'
#' @examples
#' data(IAM_argum_1984)
#' IAM.args_scenar(IAM_argum_1984)
#'
#' # Activate scenario "Sc_saisine"
#' IAM_argum_1984 <- IAM.editArgs_Scenar(IAM_argum_1984, selected = "Sc_saisine")
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
#' @rdname IAM.editArgs_Scenar
#' @export
setGeneric("IAM.editArgs_Scenar", function(object, ...){
  standardGeneric("IAM.editArgs_Scenar")
}
)

#' @name IAM.editArgs_Scenar
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
