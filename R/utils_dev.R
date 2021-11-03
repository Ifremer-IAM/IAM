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
