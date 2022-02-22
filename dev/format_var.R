#' Format IAM.output variable as a tibble.
#'
#' Extract an element from \code{\link[IAM]{iamOutput-class}} and format this as a
#' tibble. This tibble is 6 column wide with specified format (see details)
#'
#' @param name name of the exported variable from IAM. chr.
#' @param object \code{\link[IAM]{iamOutput-class}} object created with
#' \code{IAM.model()}.
#'
#' @details The reconcilSPP variable return character value,
#' so it can't be binded with other variable which are numeric.
#'
#' Format : long format for ggplot and dplyr analysis to facilitate filter and
#' summary computations. If the variable is not defined for a specific dimension
#' (column), this column is filled with NA.
#' \describe{
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
#' @name format_var
#' @rdname format_var
#'
#' @author Maxime Jaunatre
#'
#' @export
format_var <- function(name, object){

  nsp <- names(object@outputSp)
  ns <- names(object@output)


  if(name %in% ns){
    res <- format_vareco(name, object)
  } else if(name %in% nsp){
    res <- format_varsp(name, object)
  } else {
    warning(paste(name,"variable does not exist in an IAM output."),
            call. = FALSE)
    res <- NULL
  }

  return(res)
}


#' @importFrom stats ftable
#' @importFrom tidyr pivot_longer
#' @import dplyr
#' @import tibble
#' @importFrom magrittr %>%
#'
#'
#' @name format_var
#' @rdname format_var
#'
#' @export
format_vareco <- function(name, object){

  var = object@output[[name]]

  cols <- c(fleet = NA_character_, metier = NA_character_,
            age = NA_character_, year = NA_character_)

  if(!is.null(var)){
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
      add_column(., !!!cols[setdiff(names(cols), names(.))]) %>%
      .[, c("fleet", "metier", "age", "year", "value")] %>%
      add_column(variable = name, species = NA, .before = "fleet") %>%
      # mutate_if(is.factor, as.character) %>%
      mutate( value = as.numeric(value),year = as.numeric(year) )

  } else {
    res <- NULL
  }

  return(res)
}

#' @importFrom stats ftable
#' @importFrom tidyr pivot_longer
#' @import dplyr
#' @import tibble
#' @importFrom magrittr %>%
#'
#'
#' @name format_var
#' @rdname format_var
#'
#' @export
format_varsp <- function(name, object){

  var = object@outputSp[[name]]

  cols <- c(fleet = NA_character_, metier = NA_character_,
            age = NA_character_, year = NA_character_)

  len <- length(var)
  if(len == 0 | any(unlist(lapply(var, is.null)))){ return(NULL)}

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
          add_column(., !!!cols[setdiff(names(cols), names(.))]) %>%
          .[, c("fleet", "metier", "age", "year", "value")] %>%
          add_column(variable = name, species = names(var)[sp],
                     .before = "fleet")

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
      add_column(variable = name, .before = "species") %>%
      add_column(., !!!cols[setdiff(names(cols), names(.))]) %>%
      .[, c("variable", "species", "fleet", "metier", "age", "year", "value")]

  }

  res <- res  %>%
    # mutate_if(is.factor, as.character) %>%
    mutate( value = as.numeric(value),year = as.numeric(as.character(year)) )

  return(res)
}
