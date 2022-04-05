#' Format IAM.output variable as a tibble.
#'
#' Extract an element from \code{\link[IAM]{iamOutput-class}} and format this as a
#' tibble. This tibble is 6 column wide with specified format (see details)
#'
#' @param object \code{\link[IAM]{iamOutput-class}} object created with
#' \code{IAM.model()}.
#' @param name names of the exported variable from IAM. "summary" will select
#' few variables that are usually displayed for most common scenarii ("Fbar",
#' "SSB", "L_et", "N", "nbv_f", "effort2_f", "GVLav_f",
#' "gva_f", "gp_f", "wageg_f", "wagen_f"). chr.
#' @param sim_name Name of the simulation, used for joining multiple scenarii
#' iamOutput. NA by default, chr.
#' @param n Number of the simulation, used for joining replicated iamOutput.
#' NA by default, dbl.
#'
#' @details The reconcilSPP variable return character value,
#' so it can't be binded with other variable which are numeric.
#'
#' @return
#' Format : long format for ggplot and dplyr analysis to facilitate filter and
#' summary computations. If the variable is not defined for a specific dimension
#' (column), this column is filled with NA.
#' \describe{
#'   \item{sim_name}{Simulation name. chr.}
#'   \item{n}{Number of the simulation. dbl}
#'   \item{variable}{Single value repeated, but is needed if multiple variables
#'   are assembled with \code{bind_row()}. chr vector}
#'   \item{species}{Species names, can contain dynamic and static species
#'   depending on the variable selected.. chr vector}
#'   \item{fleet}{Fleet names. fct vector}
#'   \item{metier}{Metier names. fct vector}
#'   \item{age}{Ages for dynamic species. fct vector}
#'   \item{year}{year step for the simulation. dbl vector}
#'   \item{value}{Value of the variable for every column place indicated on the
#'   same row. dbl vector}
#' }
#'
#'
#' @examples
#' data("IAM_input_1984")
#' data("IAM_argum_1984")
#' sim_statu_quo <- IAM::IAM.model(objArgs = IAM_argum_1984, objInput = IAM_input_1984)
#'
#' #' cas vide return NULL
#' IAM.format(sim_statu_quo, "not a variable")
#' #' cas simple bio
#' IAM.format(sim_statu_quo, "SSB")
#' #' cas complexe bio
#' IAM.format(sim_statu_quo, "F")
#' IAM.format(sim_statu_quo, "N_S1M1")
#' #' cat eco
#' IAM.format(sim_statu_quo, name = "ratio_gp_K_f")
#' #' Warning here values are string !
#' IAM.format(sim_statu_quo, name = "reconcilSPP")
#'
#'
#'
#' @name IAM.format-method
#' @rdname IAM.format-method
#' @aliases IAM.format
#'
#' @author Maxime Jaunatre
#'
#' @export
setGeneric("IAM.format", function(object, name, sim_name = NA_character_,
                                  n = NA_real_){ # IAM.format gen ####
  standardGeneric("IAM.format")
},
signature = c("object", "name")
)

#' @rdname IAM.format-method
#' @aliases IAM.format,iamOutput-method
setMethod(
  "IAM.format", signature("iamOutput", "character"),
  function(object, name, sim_name = NA_character_, n = NA_real_){ ## iamOutput meth ####

    if(any(name == "summary")){
      name <- c(
        "Fbar", "SSB", "L_et",
        # "ratio_L_et_tac",
        "N", #BIO (N = courbe d'age)
        "nbv_f", "effort2_f", "GVLav_f", "gva_f", "gp_f", "wageg_f", "wagen_f"#ECO
      )
    }
    l <- length(name)

    if(l > 1){

      res <- vector(mode = "list", length = l)
      for(i in seq_along(name)){
        res[[i]] <- IAM.format(object, name[i], sim_name, n) # /!\ Recursive call !
      }
      res <- do.call(rbind, res)
    } else {

      ns <- names(object@output)
      nsp <- names(object@outputSp)

      if(name %in% ns){
        res <- format_vareco(object, name, sim_name, n)
      } else if(name %in% nsp){
        res <- format_varsp(object, name, sim_name, n)
      } else {
        warning(
          paste0("\"", name,"\" variable does not exist in an IAM output."),
          call. = FALSE
          )
        res <- NULL
      }
    }

    if(!is.null(res)){
      class(res) <- c("iam_formtbl", class(res))
    }

    return(res)
  }
)


#' @importFrom stats ftable
#' @importFrom tidyr pivot_longer
#' @import dplyr
#' @import tibble
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#'
#' @rdname IAM.format-method
#'
#' @export
format_vareco <- function(object, name, sim_name = NA_character_, n = NA_real_){

  var = object@output[[name]]

  cols <- c(fleet = NA_character_, metier = NA_character_,
            age = NA_character_, year = NA_character_)

  if(is.null(var)){
    return(NULL)
  }

  att <- attributes(var)$DimCst
  nam_to <- c("fleet","metier","age","year")[att>0]

  if(sum(att > 1) > 2 ){ # case with multiple dimensions

    res <- as_tibble(ftable(var))
    names(res) = c(nam_to, "value")

  } else { # double dimension

    res <- as.data.frame(var) %>%
      rownames_to_column(var = nam_to[1]) %>%
      pivot_longer( !nam_to[1], names_to = nam_to[2])
  }

  res <- res %>%
    add_column(!!!cols[setdiff(names(cols), names(res))]) %>%
    select(.data$fleet, .data$metier, .data$age, .data$year, .data$value) %>%
    add_column(sim_name = sim_name, n = n, variable = name,
               species = NA, .before = "fleet")
  # mutate_if(is.factor, as.character) %>%

  if(name == "reconcilSPP"){
    warning(paste("reconcilSPP variable has character for value.",
                  "This will cause problems if you need to rbind ouputs",
                  "of IAM.format()"))
    res <- res %>% mutate(year = as.numeric(.data$year) )
  } else {
    res <- res %>%
      mutate( value = as.numeric(.data$value),
              year = as.numeric(.data$year) )
  }



  return(res)
}

#' @importFrom stats ftable
#' @importFrom tidyr pivot_longer
#' @import dplyr
#' @import tibble
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#'
#' @rdname IAM.format-method
#'
#' @export
format_varsp <- function(object, name, sim_name = NA_character_, n = NA_real_){

  var = object@outputSp[[name]]

  cols <- c(fleet = NA_character_, metier = NA_character_,
            age = NA_character_, year = NA_character_)

  len <- length(var)
  if(len == 0 | all(unlist(lapply(var, is.null)))){ return(NULL)}

  tmp <- vector(mode = "list", length = len)
  simple <- FALSE

  for(sp in 1:len){
    x <- var[[sp]]
    if(!is.null(x)){
      att <- attributes(x)$DimCst
      if(sum(att > 1) > 1 ){ # case with multiple dimensions

        tmp[[sp]] <- as_tibble(ftable(x))
        names(tmp[[sp]]) = c(c("fleet","metier","age","year")[att>0], "value")

        tmp[[sp]] <- tmp[[sp]] %>%
          add_column(!!!cols[setdiff(names(cols), names( tmp[[sp]]))]) %>%
          select(.data$fleet, .data$metier, .data$age, .data$year, .data$value) %>%
          add_column(sim_name = sim_name, n = n, variable = name,
                     species = names(var)[sp], .before = "fleet")

      } else { # single dimension
        simple <- TRUE
        break()
      }
    }
  }

  if(!simple){ # case with multiple dimensions

    res <- do.call(rbind, tmp)

  } else { # single dimension

    att <- attributes(var[[1]])$DimCst
    res <- do.call(rbind, var) %>%
      as.data.frame() %>%
      rownames_to_column(var = "species") %>%
      pivot_longer(!"species",
                   names_to = switch(which(att > 1),
                                     "fleet","metier","age","year")) %>%
      add_column(sim_name = sim_name, n = n, variable = name,
                 .before = "species")

    res <- res %>%
      add_column( !!!cols[setdiff(names(cols), names(res) )]) %>%
      select(.data$sim_name, .data$n, .data$variable, .data$species,
        .data$fleet, .data$metier, .data$age, .data$year, .data$value)

  }

  res <- res  %>%
    # mutate_if(is.factor, as.character) %>%
    mutate( value = as.numeric(.data$value),
            year = as.numeric(as.character(.data$year)) )

  return(res)
}


#' Simplify multiple simulations
#'
#' Simplify multiple simulations and group them after computing
#' quantile and median value (per species, fleet, metier, age, year).
#' A selected simulation is kept to illustrate a possible way.
#'
#' @param var_format Tibble or data.frame format produced by
#' \code{\link[IAM]{IAM.format}}.
#' @param probs Quantile arguments for the values distribution once grouped.
#' Vector of 2 numeric values in [0, 1].
#' @param select_indiv Single numeric value which allow to conserve a single
#' simulation value for example. Default is 1.
#'
#' @return
#' Format : long format for ggplot and dplyr analysis to facilitate filter and
#' summary computations. If the variable is not defined for a specific dimension
#' (column), this column is filled with NA.
#' \describe{
#'   \item{sim_name}{Simulation name. chr.}
#'   \item{variable}{Single value repeated, but is needed if multiple variables
#'   are assembled with \code{bind_row()}. chr vector}
#'   \item{species}{Species names, can contain dynamic and static species
#'   depending on the variable selected.. chr vector}
#'   \item{fleet}{Fleet names. fct vector}
#'   \item{metier}{Metier names. fct vector}
#'   \item{age}{Ages for dynamic species. fct vector}
#'   \item{year}{year step for the simulation. dbl vector}
#'   \item{quant1}{First quantile from \code{probs} argument.}
#'   \item{quant2}{Second quantile from \code{probs} argument.}
#'   \item{median}{Median computed from the grouped values.}
#'   \item{value}{Value of the variable for every column place indicated on the
#'   same row. dbl vector}
#' }
#'
#' @examples
#' library(dplyr)
#' library(magrittr)
#'
#' data("IAM_input_1984")
#' data("IAM_argum_1984")
#' sim_statu_quo <- IAM::IAM.model(objArgs = IAM_argum_1984, objInput = IAM_input_1984)
#' res1 <- IAM.format(sim_statu_quo, c("SSB"), n = 1) %>%
#'   dplyr::filter(species == "ARC", year < 1987)
#' res2 <- mutate(res1, n = 2, value = value + rnorm(1, sd = 100))
#' res <- rbind(res1, res2)
#'
#' IAM.format_quant(res, c(0.025, 0.975), 2)
#'
#' @import dplyr
#' @import tibble
#' @importFrom magrittr %>%
#' @importFrom stats quantile median
#'
#' @author Maxime Jaunatre
#'
#' @name IAM.format_quant
#' @rdname IAM.format_quant
#'
#' @export
IAM.format_quant <- function(var_format, probs = c(.25, .75), select_indiv = 1){

  if(length(probs) > 2){
    warning("probs maximum length is 2")
  }

  if(var_format %>% select(- .data$value) %>% duplicated() %>% any()){
    stop(paste("Multiple rows have the same group attributes.",
               "This is most likely an indication of missing n information"))
  }

  if(all(is.na(var_format$n))){
    stop(paste("Simulation number must be provided for this function.",
               "Please check IAM.format() usage."))
  }
  if(! select_indiv %in% range(var_format$n)){
    warning(sprintf(
      paste("You selected a simulation number that don't exist.",
            "The maximum value is %d."),
      select_indiv <- max(var_format$n)))
  }

  if(!inherits(var_format, "iam_formtbl")){
    stop("var_format must be a table created by IAM.format().")
  }

  res <- var_format %>%
    group_by(.data$sim_name, .data$variable, .data$species,
             .data$fleet, .data$metier, .data$age, .data$year) %>%
    summarize(quant1 = quantile(.data$value, probs = probs[1], na.rm = TRUE),
              quant2 = quantile(.data$value, probs = probs[2], na.rm = TRUE),
              median = median(.data$value, na.rm = TRUE),
              sel_val = .data$value[.data$n == select_indiv],
              nsim = length(unique(n)),
              .groups = "keep") %>%
    rename(value = .data$sel_val) %>%
    ungroup()

  class(res) <- c("iam_quantbl", class(res))

  return(res)
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# Methode composite issue de IAM.format, ajoutant le regroupement des
# scenarii --> input des m?thodes graphiques

#' Formating method for grouping various \code{\link[IAM]{iamOutput-class}} or \code{\link[IAM]{iamOutputRep-class}} objects
#'
#' @param listObj A list containing \code{\link[IAM]{iamOutput-class}} OR \code{\link[IAM]{iamOutputRep-class}} objects. Mixed lists are not allowed.
#' @param ... Further arguments :
#' \describe{
#' \item{elmt}{ Name of the operating variable. chr.}
#' \item{spp}{Name of the considered species (relevant only for some variables). chr}
#' \item{agg}{Index(es) describing dimensions on which agregation must be done,
#' referring to 'DimCst' attributes (1=fleet, 2=metier, 3=age/category, 4=time). num.}
#' \item{headers}{Optionnal. Dimension that will be developped as headers in the output data.frame. chr.}
#' \item{subs}{An optionnal 'subset' argument to operate a final subsetting on output dataframe. logic}
#' }
#'
#' @rdname IAM.unite
#' @export
setGeneric("IAM.unite", function(listObj,...){  # IAM.unite gen ####
	standardGeneric("IAM.unite")
	}
)

#' @rdname IAM.unite
#' @aliases IAM.unite,list-method
setMethod("IAM.unite", signature(listObj="list"), function(listObj,...){  ## list meth ####

	if (!(all(unlist(lapply(listObj,class))%in%"iamOutput") | all(unlist(lapply(listObj,class))%in%"iamOutputRep")))
    stop("only 'iamOutput' objects OR (not AND) 'iamOutputRep' objects allowed in input list!!")



  df <- do.call("rbind",lapply(listObj, function(x) {
                                  if (class(x)%in%"iamOutput") {
                                    return(cbind.data.frame(IAM.format(x,...),
                                                            Scen=ifelse(x@arguments$Scenario$active==0,
                                                                        "Status Quo",
                                                                        x@arguments$Scenario$ALLscenario[x@arguments$Scenario$SELECTscen])))
                                  } else {
                                    ll <- IAM.format(x,...)
                                    return(cbind.data.frame(do.call("rbind",ll),
                                                            iter=rep(1:length(ll),unlist(lapply(ll,nrow))),
                                                            Scen=ifelse(x@arguments$Scenario$active==0,
                                                                        "Status Quo",
                                                                        x@arguments$Scenario$ALLscenario[x@arguments$Scenario$SELECTscen])))
                                  }

                                  }))


    call <- match.call()
    llhead <- as.list(call)[-1]$headers[1]


    if (!is.null(llhead)) {
      if (llhead%in%"Scen") {

    if ((ncol(df)>2) & (nrow(df)>0)) {
    mat <- tapply(df[,1],list(do.call('paste',c(df[,2:(ncol(df)-1),drop=FALSE],list(sep=":-:-:"))),df[,"Scen"]),function(x) x)
    mat1 <- do.call("rbind",lapply(rownames(mat),function(x) strsplit(x,":-:-:")[[1]]))
    df1 <- as.data.frame(mat1) ; names(df1) <- names(df)[-(c(1,ncol(df)))]
    df2 <- as.data.frame(mat) ; rownames(df2) <- NULL
    df <- cbind.data.frame(df1,df2)
  }

  rownames(df) <- 1:nrow(df)

    }
  }

  return(df)

})



#:::::::::::::::
#Examples
#:::::::::::::::


#arg2 <- IAM.args(arg)
#
#mod2 <- IAM.model(arg2,out)
#
#IAM.unite(list(mod,mod2),elmt = "C", spp = "Langoustine", agg = 1:2, headers = "a_c", t%in%(2012:2025))




