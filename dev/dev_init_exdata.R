
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
library(beepr)
library(openxlsx)
if(!require(IAM)) devtools::load_all()

# Objective ####
#' Modify xlsx with most of RData to obtain a clean file

tea_br <- TRUE # cancel long computations
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
  if(!dir.exists("dev/data/fleets")){
    dir.create("dev/data/fleets")
  }
  fleet <- "dev/raw_data/fleets/"
  if(dir.exists(fleet)) {
    warning(paste(fleet, "already exists and previous files won't be purged"))
  }
  dir.create(fleet, showWarnings = FALSE)
  file.copy(from = paste0(fleet,enigmaF$origin, ".csv"),
            to = paste0("dev/data/fleets/", enigmaF$code, ".csv"),
            overwrite = TRUE)
  for(i in 1:nrow(enigmaF)){
    cat(" -",enigmaF$code[i])
    nme <- paste0("dev/data/fleets/",enigmaF$code[i], ".csv")
    sheet <- read.csv(nme, sep = ";")
    sheet$annee <- 1984
    sheet <- seek_replace(sheet, enigma)
    # multiply vessel
    pattern <- grepl("nbv_f|Lref_f|H_f|rep_f|fixc_f|dep_f|ic_f|K_f|persc_f",
                     sheet$nom_variable)
    sheet$indicateur[pattern] <- as.character(
      nrep[i] * as.numeric(sheet$indicateur[pattern])
    )
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
  rowNames=FALSE,colNames=FALSE, startRow = 100, skipEmptyRows = FALSE)
  ndim <- nrow(scenar) + 99 # old place for the scenarii, to clean
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
  deleteData(wb, "Scenarii", rows = 100:ndim, cols = 1:ncol(scenar),
             gridExpand = TRUE)
  addStyle(wb, "Scenarii", createStyle(fgFill = NULL),
           rows = 100:ndim, cols = 1:ncol(scenar), gridExpand = TRUE)
  writeData(wb, sheet = "Scenarii", scenar, startRow = 100, colNames = F)
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
    writeData(wb, sheet = sh, sheet, colNames = F)
    gc()
  }
  cat("\n-done-\n")
  rm(sh, sheet, old, place, replaced_names)
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## Export dataset ####
  worksheetOrder(wb) <- match(save_names, names(wb))
  saveWorkbook(wb,"dev/data/inputFile.xlsx",overwrite = TRUE)
}
rm(wb, Mbio, Meco, nrep, save_names, enigmaSp, enigma, rawfilen, sp2rm, enigt)

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
  save(list = ls(), file = "dev/data/inpSS3darwiniana1984.RData", envir = ss3)
})
rm(ss3, enigmaF, enigmaM)




# IAM.input ####
if(tea_br){
  load("dev/data/inpSS3darwiniana1984.RData")
  input1984 <- IAM::IAM.input(fileIN = "dev/data/inputFile.xlsx",t_init=1984,nbStep=2,folderFleet="dev/data/fleets",
                              Fq_i=list(DAR=iniFq_i),iniFq_i=list(DAR=iniFq_i),Fq_fmi=list(DAR=iniFq_fmi),iniFq_fmi=list(DAR=iniFq_fmi),
                              FqLwt_i=list(DAR=iniFqLwt_i),iniFqLwt_i=list(DAR=iniFqLwt_i),FqLwt_fmi=list(DAR=iniFqLwt_fmi),iniFqLwt_fmi=list(DAR=iniFqLwt_fmi),
                              FqDwt_i=list(DAR=iniFqDwt_i),iniFqDwt_i=list(DAR=iniFqDwt_i),FqDwt_fmi=list(DAR=iniFqDwt_fmi),iniFqDwt_fmi=list(DAR=iniFqDwt_fmi),
                              Nt0s1q=list(DAR=Nt0s1q),Ni0q=list(DAR=Ni0q),iniNt0q=list(DAR=iniNt0q),matwt=list(DAR=mat_morphage),
                              verbose = TRUE)

  if(exists("input1984")) {
    beep(5) ; cat("You rock !\n ")
    save(input1984, file = "dev/data/inputIFR.RData")
  } else {
    beep(9) ; cat("Fail \n")
  }
} else {
  load("dev/data/inputIFR.RData")
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# argumIFR <- IAM.args(inputIFR)
#
# #recrutement corinna 1985 = 16402000
#
# #module Eco DCF ('active' et 'type' a priori inutiles puisque le code C++ n'int?gre plus le mod?le complet)
# argumIFR@arguments$Eco$active <- as.integer(1)
# argumIFR@arguments$Eco$type <- as.integer(2)
# argumIFR@arguments$Eco$dr <- 0.04
# argumIFR@arguments$Eco$perscCalc <- as.integer(1)
# #Gestion d?sactiv?
# argumIFR@arguments$Gestion$active <- as.integer(0)
# argumIFR@arguments$Gestion$delay <- as.integer(1)
# argumIFR@arguments$Gestion$mfm[] <- (inputIFR@input$Fleet$effort1_f_m*inputIFR@input$Fleet$effort2_f_m*inputIFR@input$Fleet$nbv_f_m)/
#   as.vector(inputIFR@input$Fleet$effort1_f*inputIFR@input$Fleet$effort2_f*inputIFR@input$Fleet$nbv_f)
# argumIFR@arguments$Gestion$mfm[is.na(argumIFR@arguments$Gestion$mfm)] <- 0
# #Sc?nario d?sactiv?
# argumIFR@arguments$Scenario$active <- as.integer(0)
#
# save(x=argumIFR,file='Z:/Projet/PG GG/These_Florence/Parametrage_IAM/3.PARAMETRAGE/argumIFR_newStocks_MNZstatic_20yr.RData')

