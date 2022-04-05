#' Script for development procedure.
#'
#' Developping IAM as a package depend on {devtools} package.
#' It also include gitlab or github CI and {pkgdown} package
#' for documentation presentation.
#'
#' Below are some of the procedure to run once an edit has been
#' done in the package.
#'
#' Most of them need to be used before git commit.
#'
#' Run this script interactively and only select some chunks !
#'
#' maxime Jaunatre 31/03/2022

library(usethis)
library(devtools)
library(pkgdown)
library(fs)

# Command for pkg compil ####

## Reload IAM in empty env. ####
devtools::unload("IAM")
.rs.restartR()
devtools::load_all()

## Force C++ compil without edition ####
devtools::unload("IAM")
.rs.restartR()
devtools::clean_dll()
devtools::load_all()

## Reinstall package ####
devtools::document()
devtools::unload("IAM")
.rs.restartR()
remove.packages("IAM")
devtools::install(upgrade = "never",
                  build_vignettes = TRUE)

## Check the full package ####
devtools::check() # takes few minutes

## Compile documentation site
devtools::load_all()
rmarkdown::render("README.rmd", clean = TRUE, quiet = TRUE)
fs::file_delete("README.html")
pkgdown::build_site(preview = TRUE) # will open website

## Clear objects in env ####
rm(list = ls())
gc()
# Restart R
.rs.restartR()

#' I tested using some options but it's poorly used in the package.
# options(dev = TRUE)
# options(IAM.dev = TRUE) # TRUE show verbose, FALSE not. NULL equal FALSE

# Index of files ####
#' @name dev_cp_datarmor File copy and name removing from datarmor runs
#' in F. work (2017-2018). Require access to Z: disk. Copy files in
#' dev/raw_data/. File exported in gitlab IAM/dataex_setup.
#'
#' @name dev_graphic Work file to test plots about IAM.
#' Objective is to remove it and present plots in main vignette.
#'
#' @name dev_init_exdata Script to build IAM_xxx_1984 example dataset.
#' It use files in dev/raw_data/ and produce cache dev/data .RData to be
#' more efficient.
#'
#' @name dev_init_simpl_exdata Similar to previous file but aim to build a
#' single .xlsx file (rm fleet dir dependancy and SS3 species.). Not finalised.
#'
#' @name dev_load_oldRdata Old .RData produced by IAM are dependent of these
#' versions because S4 classes remember package of creation.
#' Tricks R but creating an empty namespace. To use with caution.
#'
#' @name format_var Scrip to test IAM.model, IAM.format and IAM.format_quant
#' @name shiny_plot Scrip for testing the IAM.test_plot function
#' @name shiny_test Script for testing the IAM.args function and shiny app.
#'
#' @name dev_manhattan Small exploration project for IAM.input refactoring
#' and speed
#' @name stripe_input Small exploration project for IAM.input refactoring
#' and speed

# Notes ####
#' I can't put simulation in data/ because maximum recommended is 5MB
