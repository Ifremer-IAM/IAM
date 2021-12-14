

#' @importFrom tidyr pivot_longer pivot_wider
#' @param file filepath with an .xlsx extention.
#' @param sheet sheet name.
# TODO : is this very usefull since we always load specific one ?
read_fm_mat <- function(file, sheet){
  prefix <- c("v__","t__","i__","f__","m__","l__","e__","c__")
  FM <- read.xlsx(file, sheet = sheet, rowNames=FALSE, colNames=FALSE,
                  skipEmptyRows = FALSE, skipEmptyCols = FALSE)
  FM[is.na(FM)] <- ""
  # cut table in list
  indicRow <- apply(FM,1,function(x) any(substring(x,1,3) %in% prefix))
  indicTbl <- cumsum(apply(FM,1,function(x) all(x=="")))
  sepTabl <- split(FM[indicRow,],indicTbl[indicRow])
  # format and merge tables
  tables <- lapply(sepTabl, function(x){
    var <- as.character(x[1,grepl("m__", x[1,])])
    nt <- x[-1,] ; nt <- nt[, !apply(nt,2, function(x) all(x ==""))]
    colnames(nt) <- c("X", "species", "fleet", var)
    nt <- pivot_longer(nt, cols = -c(1:3))
    return(nt)
  })
  # require to merge all and pivot the other way and filter
  tables <- pivot_wider(do.call(rbind, tables), values_fill = "")

  tables <- split(tables, tables$species)


  tables <- lapply(tables, function(x){
    fid <- sub("^f__", "", x$fleet)
    Mbio <- sub("^m__", "", colnames(x))[-c(1:3)]
    x <- apply(x[,-c(1:3)], 2, as.numeric)
    dimnames(x) <- list(fid, Mbio)
    attr(x, "Dimcst") <- c(dim(x), 0, 0)
    return(x)
  })
  return(tables)
}

read_icat_mat <- function(file, sheet){
  prefix <- c("v__","t__","i__","f__","m__","l__","e__","c__")
  icat <- read.xlsx(file, sheet = sheet, rowNames=FALSE, colNames=FALSE,
                    skipEmptyRows = FALSE, skipEmptyCols = FALSE)
  icat[is.na(icat)] <- ""
  # cut table in list
  indicRow <- apply(icat,1,function(x) any(substring(x,1,3) %in% prefix))
  indicTbl <- cumsum(apply(icat,1,function(x) all(x=="")))
  sepTabl <- split(icat[indicRow,],indicTbl[indicRow])
  # format and merge tables
  tables <- lapply(sepTabl, function(x){
    var <- as.character(x[1,grepl("c__", x[1,])])
    nt <- x[-1,-1] ; nt <- nt[, !apply(nt,2, function(x) all(x ==""))]
    colnames(nt) <- c("species", "age", var)
    return(nt)
  })
  names(tables) <- lapply(tables, function(x) unique(x$species))

  tables <- lapply(tables, function(x){
    age <- sub("i__", "", x$age)
    cat <- sub("c__", "", colnames(x))[-c(1:2)]
    x <- apply(x[,-c(1:2)], 2, as.numeric)
    dimnames(x) <- list(age, cat)
    return(x)
  })
  return(tables)
}

read_mm_mat <- function(file, sheet){
  MM <- read.xlsx(file, sheet = "mm_matrix", rowNames=FALSE,
                  skipEmptyRows = FALSE, skipEmptyCols = FALSE)[,-1]
  colnames(MM) <- c(sub("espèce", "species", colnames(MM)))
  MM <- split(MM, MM$species)

  MM <- lapply(MM, function(x){
    fid <- paste0(sub("^f__","",x$flottille), x$'mStock\\mFleet')
    Mbio <- sub("^m__", "", colnames(x)[-c(1:3)])
    x <- apply(x[,-c(1:3)], 2, as.numeric)
    dimnames(x) <- list(fid, Mbio)
    return(x)
  })
  return(MM)
}

#' @param Market Market table
#' @param Sp Species label with format "e__XXX"
get_cat <- function(Market, Sp){
  Sp_cat  <- unique(Market$categorie[Market$espece == Sp])
  Sp_cat <- Sp_cat[Sp_cat != ""]
}


#' mock function to refactor the read.input
#' @importFrom openxlsx getSheetNames read.xlsx
#' @importFrom methods rbind2 new
#' @importFrom utils read.table
stripe.input <- function(file, t_init, nbStep, t_hist_max = t_init,
                       desc = "My input", folderFleet = NULL, verbose = FALSE ) {
  rm(list = ls())
  library(openxlsx)
  file <- "inst/extdata/IAM_MED_simpl.xlsx"
  t_init <- 2020
  t_hist_max <- 2020
  nbStep <-  5
  desc <- "Med tryhard"
  folderFleet = NULL
  verbose = TRUE

  ## Init return ####
  specific <- list(Species=NULL, StaticSpp=NULL, AllSpp = NULL,
                   Fleet=NULL, Metier=NULL, MetierEco=NULL,
                   Ages=NULL, Cat=NULL,
                   t_init=t_init, NbSteps=as.integer(nbStep),
                   times=as.integer(seq(t_init,length=nbStep)),
                   Q=NULL, S=NULL)
  historical <- list()
  input <- list()
  scenario <- list()
  stochastic <- list()

  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::READ_ALL
  ## Read sheets ####
  # Also build specific !
  tbls <- getSheetNames(file)
  prefix <- c("v__","t__","i__","f__","m__","l__","e__","c__")
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::Stocks
  if(verbose) cat("Read stocks\n") ### Stocks ####
  namSp <- tbls[grepl("Stock__", tbls)]
  nSp <- length(namSp)
  if(nSp == 0) stop("Missing Dynamic stock !")
  # extract names to specific
  cl_namSp <- sub("Stock__", "", namSp)
  AgeCat <- vector(mode = "list", length = nSp)
  SQ <- integer(length = nSp)
  specific$Species <- names(AgeCat) <- names(SQ) <- cl_namSp
  specific$Ages <- specific$Cat <- AgeCat
  specific$Q <- specific$S <- SQ
  rm(SQ, AgeCat, cl_namSp)

  # require tbls, file, namSp
  Stock <- vector(mode = "list", length = nSp)
  names(Stock) <- namSp
  for(stock in namSp){
    Stock[[stock]] <- read.xlsx(file, sheet = stock, rowNames=FALSE, colNames=FALSE,
                            skipEmptyRows = FALSE, skipEmptyCols = FALSE)
    Stock[[stock]][is.na(Stock[[stock]])] <- ""

    recode <- Stock[[stock]][5:34,1:4]
  }
  # return Stock
  tbls <- tbls[!grepl("^Stock__", tbls)]


  #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::Fleet
  if(verbose) cat('Read fleet : ') ### Fleets ####
  # require folderFleet, tbls, verbose
  if(!is.null(folderFleet)){ # Force on fleet folder
    namF <- sort(list.files(folderFleet, pattern = "^f__.*.csv$"))
    if(length(namF) == 0) stop("Missing fleet !")
    Fleet <- vector(mode = "list", length = length(namF))
    names(Fleet) <- sub(".csv$", "", namF)
    for(f in namF){
      Fleet[[sub(".csv$", "", f)]] <- read.csv(file.path(folderFleet, f),
                                               sep = ";")
      if(verbose) cat(' .')
    }
  } else {
    namF <- sort(tbls[grepl("^f__", tbls)])
    if(length(namF) == 0) stop("Missing fleet !")
    Fleet <- vector(mode = "list", length = length(namF))
    names(Fleet) <- namF
    for(f in namF){
      Fleet[[f]] <- read.xlsx(file, sheet = f, rowNames=FALSE,
                              skipEmptyRows = FALSE, skipEmptyCols = FALSE)
      # replace NA with "" for chr cols.
      Fleet[[f]] <- as.data.frame(sapply(Fleet[[f]], function(x) {
        if(is.character(x)){ x[is.na(x)] <- "" } ; return(x)
      }, simplify = FALSE))
      if(verbose) cat(' .')
    }
  }
  rm(f, folderFleet) # return namF and Fleet list
  ### Extract MetierEco and AllSpp ####
  specific$Fleet <- sub("f__", "", namF)
  MetierEco <- unlist(lapply(Fleet, function(x) unique(x$metier)))
  MetierEco <- unique(MetierEco[MetierEco != ""])
  specific$MetierEco <- sub("m__", "", MetierEco)
  AllSpp <- unlist(lapply(Fleet, function(x) unique(x$espece)))
  AllSpp <- unique(AllSpp[AllSpp != ""])
  specific$AllSpp <- sub("e__", "", AllSpp)
  specific$StaticSpp <- specific$AllSpp[! specific$AllSpp %in% specific$Species]


  tbls <- tbls[!grepl("^f__", tbls)]
  if(verbose) cat(' OK\n')


  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::Market
  if(verbose) cat('Read Market') ### Market ####
  # require verbose, file
  Market <- read.xlsx(file, sheet = "Market", rowNames=FALSE,
                      skipEmptyRows = FALSE, skipEmptyCols = FALSE)
  Market <- as.data.frame(sapply(Market, function(x) {
    if(is.character(x)){ x[is.na(x)] <- "" } ; return(x)
  }, simplify = FALSE))


  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::Matrix
  if(verbose) cat('Read matrix : ') ### Matrix ####
  # require file
  if(any(!c("fm_matrix", "icat_matrix", "mm_matrix") %in% tbls)) {
    stop("Missing conversion matrix (mm, fm or icat)")
  }
  FM <- read_fm_mat(file, "fm_matrix")
  ICAT <- read_icat_mat(file, "icat_matrix")
  MM <- read_mm_mat(file, "mm_matrix")
  tbls <- tbls[!grepl("_matrix", tbls)]
  if(verbose) cat(' OK\n')


  #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::Reste
  tablist <- vector(mode = "list", length = length(tbls))
  names(tablist) <- tbls
  if(verbose) cat("Reading file : ") ### Rest ####
  for(tab in tbls){
    if(verbose) cat(tab, "")
    tablist[[tab]] <- read.xlsx(file, sheet = tab, rowNames=FALSE, colNames=FALSE,
                              skipEmptyRows = FALSE, skipEmptyCols = FALSE)
    # replace NA with "" for chr cols.
    tablist[[tab]][is.na(tablist[[tab]])] <- ""
  }
  rm(tab, tbls, file, nbStep, t_hist_max, t_init)

  if(verbose) cat("OK \n")
  #:::: EOF READ sheets : tablist, fleet


  #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::INPUT
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::HISTORICAL
  ## Build historical and input with loop ####
  # names(input) <- names(historical)
  # TODO : put all matrix in input$stock


  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::SCENARIO
  ## Read and build Scenario ####

  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::STOCHASTIC
  ## Build stochastic ####

  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::OPTIMISATION
  ## Build optimisation ####




  return(new("iamInput", desc=desc, specific=specific,
             historical = historical, input = input,
             scenario = scenario, stochastic = stochastic))
}
