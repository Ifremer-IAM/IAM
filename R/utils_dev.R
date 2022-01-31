#' Impact Assesment Model
#'
#' IAM (Impact Assessment Model) is a bio-economic model
#' for the simulation of fisheries dynamics,
#' integrating specific decision support tools within the
#' framework of the theoretical implementation of management measures
#' It is a discrete time, multi-fleet, multi-trade, multi-species annual
#' bio-economic model with "age" components for the biological part,
#' and "commercial category" components for the economic part.
#'
#' @docType package
#' @name IAM
NULL


#' If x is `NULL`, return y, otherwise return x
#'
#' from ThinkR style
#'
#' @param x,y Two elements to test, one potentially `NULL`
#'
#' @noRd
#'
#'
#' @examples
#' NULL %||% 1
"%||%" <- function(x, y){
  if (is.null(x)) {
    y
  } else {
    x
  }
}

#' return `TRUE` if in `devloppment mode`
#'
#' from ThinkR style
#'
#' @export
app_dev <- function(){
  getOption( "IAM.dev" ) %||% FALSE
}

#' Initiate a project with files and directory
#'
#' This function can be used to create a good practice environment for further
#' usage of the IAM package.
#'
#' @details
#' R directory is intended to store R scripts with functions.
#' Use \code{source("R/script.R")} to load functions inside a script.
#'
#' output directory is intended to store all output created from this project.
#' All output can can be deleted at any point and obtained from the main script.
#'
#' data contains two directory for raw and processed dataset.
#' Raw data should not be modified during the code evaluation.
#' Processed data can be deleted at any point and obtained from the main script.
#' If a git repository is created before calling this function,
#' data/raw/ directory will by default be added to .gitignore file.
#'
#' @param project Title of the project in which IAM is used
#' @param file Name of the analysis file initialised with the project.
#' Used as the main script of the project.
#' NULL will cancel the file creation.
#'
#' @import fs
#' @importFrom usethis use_git_ignore
#' @importFrom rstudioapi isAvailable navigateToFile
#'
#' @author Maxime Jaunatre
#'
#' @export
IAM.dev <- function(project = "Workgroup", file = "analysis.R"){

  # cree les dossiers
  dir_to_create <- c("R", "output", "data/raw", "data/processed")

  for(dir in dir_to_create){
    if(!dir_exists(dir)){
      cat("creating", dir, "\n")
      dir_create(dir)
    } else {
      cat(dir, "already exist, check what's inside before working with it.\n")
    }
  }
  # edit .gitignore
  if(dir_exists(".git")){
    usethis::use_git_ignore(
      ignores = c(".Rproj.user", ".Rhistory", ".RData", ".Ruserdata", "data/raw/"),
      directory = "./")
  }


  # edit analyses.R
  if(!is.null(file)){
    is_proj <- any(grepl(".*[[:punct:]]Rproj$", dir_ls()))

    if(file_exists(file)){
      warnings(c(file, "already exist"))
      choice <- menu(c("Yes", "No"),
                     title = sprintf("Replace %s as default analysis script ?", file))
      if(choice == 2){
        file <- readline("New default analysis script :")
        if(!grepl(".*[[:punct:]]R$", file)){file <- paste0(file,".R")}
      }
    }
    file_create(file)

    write(c("###########################",
            sprintf("\n# IAM analysis for %s #\n", project),
            "###########################",
            sprintf("# date : %s", format(Sys.time(), "%d %b %Y"))),
          file = file)

    if(is_proj){
      write(sprintf('# Working with local files\nsetwd("%s")', getwd()),
            file = file)
    }
    write(c("\n# Loading library####\nlibrary(IAM)\n####",
            "\n# Import data####\nIAM.input()\n####", # TODO : edit this with default dataset
            "\n# Create arguments####\nIAM.args()\n####",
            "\n# Run model####\nIAM.model()\n####",
            "\n# Exploit data output####\n \n \n####"),
          file = file, append = TRUE)

    if(rstudioapi::isAvailable()) {rstudioapi::navigateToFile(file)}
  }
}


#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column, add_column, as_tibble
#'
#' @author Maxime Jaunatre
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
      add_column(variable = name, .before = "fleet")

  } else {
    res <- NULL
  }

  res <- res  %>%
    mutate( value = as.numeric(value),year = as.numeric(year) )
  return(res)
}

#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column, add_column, as_tibble
#'
#' @author Maxime Jaunatre
#'
#' @export
format_varsp <- function(name, object){

  var = object@outputSp[[name]]

  cols <- c(fleet = NA_character_, metier = NA_character_,
            age = NA_character_, year = NA_character_)

  len <- length(var)
  if(len == 0){ return(NULL)}

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
    mutate( value = as.numeric(value),year = as.numeric(year) )
  return(res)
}

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
