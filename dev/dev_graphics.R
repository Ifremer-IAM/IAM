

devtools::load_all()

load("dev/data/inputIFR.RData")
load("dev/data/argumIFR.RData")

sim1984 <- IAM::IAM.model(objArgs = argum1984, objInput = input1984, verbose = TRUE)


#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column, add_column, as_tibble
format_variable_sp <- function(name, object = sim1984,
                            specific = sim1984@specific){

  var = object@outputSp[[name]]

  cols <- c(fleet = NA_character_, metier = NA_character_,
            age = NA_character_, year = NA_character_)

  if(length(var) == 0){ return(NULL)}
  len <- length(var)
  tmp <- vector(mode = "list", length = len)
  simple <- FALSE

  for(sp in 1:len){
    x <- var[[sp]]
    if(!is.null(x)){
      att <- attributes(x)$DimCst
      if(sum(att > 1) > 1 ){ # case with multiple dimensions

        tmp[[sp]] <- as.tibble(ftable(x))
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

  return(res)
}


# cas vide PQuot
format_variable_sp("PQuot")
# cas simple SSB
format_variable_sp("SSB")
# cas complexe F
format_variable_sp("F") # %>% filter(is.na(value)) %>% head # TODO why NA here ?
format_variable_sp("N_S1M1")


