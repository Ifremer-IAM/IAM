
#=====================================================#
## How to gather old data and create the example one ##
#=====================================================#

#' Example dataset is from datarmor run from Florence (2017-2018).
#' original script is Misc.r in origin (file path defined below)
#'
#' all raw_data files are copied with dev_cp_datarmor.R.
#'
#' @author Maxime Jaunatre
#' @details 14/12/2021
#'
#' It's a xlsx file and a directory of individual fleet.
#' there is some input in .RData,
#'    * 1 for ss3
#'    * 4 icat icat matrix
#'    * 4 mortality mat
rm(list = ls())
library(beepr) ; library(praise)
library(openxlsx)
if(!require(IAM)) devtools::load_all()

# Objective ####
#' Modify xlsx with most of RData to obtain a clean file

tea_br <- TRUE # cancel long computations by setting FALSE
# Theses files are not public.
rawfile <- "dev/raw_data/inputFile.xlsx"
enigmaF <- read.csv("dev/raw_data/enigmaF.csv")
enigmaM <- read.csv("dev/raw_data/enigmaM.csv")
enigmaSp <- read.csv("dev/raw_data/enigmaSp.csv")
enigt <- enigma <- read.csv("dev/raw_data/enigma.csv")

sp2rm <- enigmaSp$origin[enigmaSp$code=="GAR"]

#' @param sheet an input.xlsx sheet for IAM
#' @param enigma dataframe with origin and code column, all character.
#' enigma tables are save in csv files in dev/raw_data
#'
#' @details ENIGMA tables a for anonymisation of dataset and should not
#' be public !
#'
#' Function that lapply in character dataframe to sub origin string with
#' code string from an enigma df.
#'
#' @author Maxime Jaunatre.
#'
seek_replace <- function(sheet, enigma){
  for(i in 1:nrow(enigma)){
    sheet <- lapply(sheet, function(x){
      sub(enigma$origin[i], enigma$code[i], x)
    })
  }
  return(as.data.frame(sheet))
}

# Simplify dataset ####
if(tea_br){
  #' Running gag for species is Coenonympha species.
  enigma <- rbind(enigmaSp, enigmaF, enigmaM, enigma)
  nrep <- c(5, 36, 24, 18, 9, 15, 60) # values by hand
  ## Modify input xlsx ####
  wb <- loadWorkbook(rawfile)
  save_names <- names(wb)
  #' selectionner les especes d'interet
  removeWorksheet(wb, paste0("Stock__", sp2rm))
  save_names <- save_names[!grepl(sp2rm, save_names)]
  replaced_names <- save_names
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## Select and rename fleet ####
  cat("-Fleets \n ")
  if(!dir.exists("inst/extdata/fleets")){
    dir.create("inst/extdata/fleets")
  } else {
    warning("dev/data/fleets already exists and previous files won't be purged")
  }
  for(i in 1:nrow(enigmaF)){
    cat(" -",enigmaF$code[i])
    nme <- paste0("dev/raw_data/fleets/",enigmaF$origin[i], ".csv")
    sheet <- read.csv(nme, sep = ";")
    sheet$annee <- 2009
    sheet <- seek_replace(sheet, enigma)
    # multiply vessel
    pattern <- grepl("nbv_f|[^V]Lref_f|H_f|rep_f|fixc_f|dep_f|ic_f|K_f|persc_f",
                     sheet$nom_variable)
    sheet$indicateur[pattern] <- as.character(
      nrep[i] * as.numeric(sheet$indicateur[pattern])
    )
    nme <- paste0("inst/extdata/fleets/", enigmaF$code[i], ".csv")
    write.csv2(sheet, nme, row.names = FALSE)
  }
  rm(i, sheet, nme, pattern)
  cat("\n-done-\n\n")
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## ICAT ####
  cat("-ICAT_matrix \n ")
  # modify the icat matrix from Florence Validation !
  ls_icat <- list.files("dev/raw_data/")[grep("newKeyini",
                                              list.files("dev/raw_data/"))]
  ls_icat <- ls_icat[!grepl(sp2rm, ls_icat)] # rm the specie sp2rm
  icats <- vector("list", length(ls_icat))
  for(i in seq(ls_icat)){
    tmp_env <- new.env()
    load(paste0("dev/raw_data/",ls_icat[i]), envir = tmp_env)
    cat(" - ", ls(tmp_env)[1])
    sp <- toupper(sub("^(newKeyini_)(.*)(_prixIndiv.RData)$", "\\2", ls_icat[i]))
    tmp <- tmp_env[[ls(tmp_env)[1]]]
    icats[[i]] <- cbind("v__icat", paste0("e__", sp), paste0("i__",rownames(tmp)),
                        tmp)
    icats[[i]] <- rbind(
      "","", "", c("", "", "", paste0("c__", colnames(tmp))),icats[[i]]
    )
  }
  rm(sp, i, tmp_env, tmp, ls_icat)

  ncols <- max(unlist(lapply(icats, ncol)))
  for(i in length(icats)){
    n <- ncol(icats[[i]])
    if(n < ncols){
      icats[[i]] <- cbind(icats[[i]], matrix(
        "", ncol = ncols - n, nrow <- nrow(icats[[i]]))
      )
    }
  }
  rm(n, i, ncols)
  icats <- do.call(rbind, icats)
  icats <- seek_replace(as.data.frame(icats), enigma)
  removeWorksheet(wb, "icat_matrix") ; addWorksheet(wb, sheet = "icat_matrix")
  writeData(wb, sheet = "icat_matrix", icats,
            colNames = FALSE, rowNames = FALSE)
  replaced_names <- replaced_names[replaced_names != "icat_matrix"]
  rm(icats)
  cat("\n-done-\n\n")
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## Market ####
  cat("-Market ")
  #' edit market sheet to rm fleets
  market <- read.xlsx(rawfile, sheet = "Market")
  market <- market[is.na(market$flottille) |
                     market$flottille %in% enigmaF$origin,]
  market <- market[!(market$espece == paste0("e__", sp2rm) &
                       !is.na(market$categorie)),]
  market <- seek_replace(market, enigma)
  removeWorksheet(wb, "Market") ; addWorksheet(wb, sheet = "Market")
  writeData(wb, sheet = "Market", market, colNames = F)
  replaced_names <- replaced_names[replaced_names != "Market"]
  rm(market, fleet)
  cat("-done-\n")
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## FM_matrix ####
  cat("-fm_matrix ")
  #' edit fm_sheet to rm fleets & filter unset metier
  fm_mat <- read.xlsx(rawfile, sheet = "fm_matrix",
                      rowNames=FALSE,colNames=FALSE, skipEmptyRows = FALSE)
  fm_mat <- fm_mat[
    (is.na(fm_mat$X4) & is.na(fm_mat$X3) ) |
    (is.na(fm_mat$X2) & fm_mat$X4 != paste0("m__", sp2rm)) | (
      fm_mat$X2 != paste0("e__", sp2rm) & fm_mat$X3 %in% enigmaF$origin),
    ]
  fm_mat <- fm_mat[
    ,
    colSums(apply(fm_mat, 2, function(x) grepl("m__", x)) |
              is.na(fm_mat)) < nrow(fm_mat)
  ]
  Mbio <- unique(unlist(fm_mat)[grep("m__", unlist(fm_mat))])
  fm_mat <- seek_replace(fm_mat, enigma)
  removeWorksheet(wb, "fm_matrix") ; addWorksheet(wb, sheet = "fm_matrix")
  writeData(wb, sheet = "fm_matrix", fm_mat, colNames = F)
  replaced_names <- replaced_names[replaced_names != "fm_matrix"]
  rm(fm_mat)
  cat("-done-\n")
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## MM_matrix ####
  cat("-mm_matrix ")
  #' edit mm_sheet to rm fleets & sp
  mm_mat <- read.xlsx(rawfile, sheet = "mm_matrix",
                      rowNames=FALSE,colNames=FALSE)
  mm_mat <- mm_mat[mm_mat$X2 %in% c("flottille", enigmaF$origin) &
                     mm_mat$X3 != paste0("e__", sp2rm),]
  mm_mat <- mm_mat[
    ,
    colSums(apply(mm_mat, 2, function(x) grepl("m__", x)) |
              mm_mat == 0, na.rm = TRUE) < nrow(mm_mat)
  ]
  Meco <- unique(unlist(mm_mat)[grep("m__", unlist(mm_mat))])
  Meco <- Meco[! Meco %in% Mbio]
  mm_mat <- seek_replace(mm_mat, enigma)
  removeWorksheet(wb, "mm_matrix") ; addWorksheet(wb, sheet = "mm_matrix")
  writeData(wb, sheet = "mm_matrix", mm_mat, colNames = F)
  replaced_names <- replaced_names[replaced_names != "mm_matrix"]
  rm(mm_mat)
  cat("-done-\n")
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## Scenarii ####
  cat("-Scenarri ")
  #' edit scenarri sheet to rm fleets & sp
  #' TODO : scenarii header should be removed, outdated !
  scenar <- read.xlsx(rawfile, sheet = "Scenarii",
  rowNames=FALSE,colNames=FALSE, startRow = 101, skipEmptyRows = FALSE)

  ndim <- nrow(scenar) + 100 # old place for the scenarii, to clean
  scenar <- scenar[
    is.na(scenar$X3) | grepl("e__", scenar$X3) | scenar$X3 %in% enigmaF$origin
    ,
  ]
  # need to remove 5 lines because no longer tables.
  # adding original t__year in col5 plant the scenarii table btw
  r <- which(!is.na(scenar$X5) & scenar$X5 == enigt$origin[1])
  scenar <- scenar[-r[is.na(scenar$X5[r+1]) & is.na(scenar$X4[r+1])], ]
  rm(r)
  scenar <- seek_replace(scenar, enigma)
  deleteData(wb, "Scenarii", rows = 101:ndim, cols = 1:ncol(scenar),
             gridExpand = TRUE)
  addStyle(wb, "Scenarii", createStyle(fgFill = NULL),
           rows = 101:ndim, cols = 1:ncol(scenar), gridExpand = TRUE)
  # writeData(wb, sheet = "Scenarii", scenar, startRow = 100, colNames = F)
  replaced_names <- replaced_names[replaced_names != "Scenarii"]
  rm(scenar, ndim)
  cat("-done-\n")
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## Rename rest ####
  # Change names in all xlsx sheets.
  cat("Replacing name in all sheets :\n")
  for(sh in replaced_names){
    cat(" -", sh)
    sheet <- read.xlsx(rawfile, sheet = sh, rowNames=FALSE,colNames=FALSE,
                       skipEmptyRows = FALSE, skipEmptyCols = FALSE)
    sheet <- seek_replace(sheet, enigma)
    place <- grep(substring(sh, 8), enigma$origin)
    if(grepl("Stock__", sh) & length(place) > 0){
      old <- sh
      sh <- sub(enigma$origin[place], enigma$code[place], sh)
      save_names[save_names == old] <- names(wb)[names(wb) == old] <- sh
    }
    removeWorksheet(wb, sh) ; addWorksheet(wb, sheet = sh)
    writeData(wb, sheet = sh, sheet, colNames = F)
    gc()
  }
  cat("\n-done-\n")
  rm(sh, sheet, old, place, replaced_names)
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## Export dataset ####
  worksheetOrder(wb) <- match(save_names, names(wb))
  saveWorkbook(wb,"inst/extdata/inputFile.xlsx",overwrite = TRUE)

  rm(wb, Mbio, Meco, nrep, save_names)
}
rm(enigmaSp, enigma, rawfile, sp2rm, enigt)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Edit SS3 input ####
ss3 <- new.env()
with(ss3, {
  load("dev/raw_data/inpSS3darwiniana1984.RData")
  var <- ls() ; var <- var[grep("_fmi$", var)]
  shfile <- apply(enigmaF, 2, substring, 4)
  mfile <- apply(enigmaM, 2, substring, 4)
  for(i in seq_along(var)){
    eval(parse(text = paste0(var[i], " = ", var[i],"[,,shfile[,'origin'],,]")))
    eval(parse(text = paste0("dimnames(", var[i], ")[[3]] = shfile[,'code']")))
    eval(parse(text = paste0("dimnames(", var[i], ")[[4]] = mfile[,'code']")))
  }
  rm(i, var, shfile, mfile)
  save(list = ls(), file = "inst/extdata/IAM_SS3_2009.RData", envir = ss3)
})


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# IAM.input ####
if(tea_br){
  load(IAM_example("IAM_SS3_2009.RData"))
  IAM_input_2009 <- IAM::IAM.input(
    fileIN = IAM_example("inputFile.xlsx"),
    t_init = 2009, nbStep = 12, folderFleet = IAM_example("fleets"),
    Fq_i = list(DAR = iniFq_i), iniFq_i = list(DAR = iniFq_i),
    Fq_fmi = list(DAR = iniFq_fmi), iniFq_fmi = list(DAR = iniFq_fmi),
    FqLwt_i = list(DAR = iniFqLwt_i), iniFqLwt_i = list(DAR = iniFqLwt_i),
    FqLwt_fmi = list(DAR = iniFqLwt_fmi), iniFqLwt_fmi = list(DAR = iniFqLwt_fmi),
    FqDwt_i = list(DAR = iniFqDwt_i), iniFqDwt_i = list(DAR = iniFqDwt_i),
    FqDwt_fmi = list(DAR = iniFqDwt_fmi), iniFqDwt_fmi = list(DAR = iniFqDwt_fmi),
    Nt0s1q = list(DAR = Nt0s1q), Ni0q = list(DAR = Ni0q),
    iniNt0q = list(DAR = iniNt0q), matwt = list(DAR = mat_morphage),
    verbose = TRUE
  )
  # Clean
  rm("Fq_fmi", "Fq_i", "FqDwt_fmi", "FqDwt_i", "FqLwt_fmi", "FqLwt_i",
       "iniFq_fmi", "iniFq_i", "iniFqDwt_fmi", "iniFqDwt_i", "iniFqLwt_fmi",
       "iniFqLwt_i", "iniNt0q", "mat_morphage", "Ni0q", "Nt0s1q")
  # Export
  if(exists("IAM_input_2009")) {
    beep(5) ; cat("\U0001f947",praise(),"\U0001f947\n")
    # save(IAM_input_2009, file = "dev/data/IAM_input_2009.RData")
    usethis::use_data(IAM_input_2009, overwrite = TRUE)
  } else {
    beep(9) ; cat("\U0001f624", "Keep trying!", "\U0001f624")
  }
} else {
  data("IAM_input_2009")
}

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# IAM.arg ####
if(tea_br){
  IAM_argum_2009 <- IAM.input2args(IAM_input_2009)
  # tmp <- IAM::IAM.args(input2009)
  # #recrutement corinna 1985 = 16402000
  #
  # #module Eco DCF ('active' et 'type' a priori inutiles puisque le code C++ n'integre plus le module complet)
  # IAM_argum_2009@arguments$Eco$active <- as.integer(1)
  # IAM_argum_2009@arguments$Eco$type <- as.integer(2)

  # IAM_argum_2009 <- IAM.editArgs_Eco(IAM_argum_2009, dr = 0.04, perscCalc = 1)

  # #Gestion desactive
  # IAM_argum_2009@arguments$Gestion$active <- as.integer(0)
  # IAM_argum_2009@arguments$Gestion$delay <- as.integer(1)
  # IAM_argum_2009@arguments$Gestion$mfm[] <- with(input2009@input$Fleet,{
    # (effort1_f_m * effort2_f_m * nbv_f_m) / as.vector(effort1_f * effort2_f * nbv_f)
  # })
  # IAM_argum_2009@arguments$Gestion$mfm[is.na(IAM_argum_2009@arguments$Gestion$mfm)] <- 0
  # #Scenario desactive
  # IAM_argum_2009@arguments$Scenario$active <- as.integer(0)
  #
  # save(IAM_argum_2009, file = "dev/data/argumIFR.RData")
  usethis::use_data(IAM_argum_2009, overwrite = TRUE)
} else {
  data("IAM_argum_2009")
}

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# IAM.model ####
devtools::load_all()
data("IAM_input_2009")
data("IAM_argum_2009")
# sim2009 <- IAM::IAM.model(objArgs = IAM_argum_2009, objInput = IAM_input_2009, verbose = TRUE)


