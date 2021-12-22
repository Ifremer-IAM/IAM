rm(list = ls())
file <- "dev/data/inputFile.xlsx"
# t_init <- 2020
# t_hist_max <- 2020
# nbStep <-  5
# desc <- "Med tryhard"
# folderFleet = "dev/data/fleets"
verbose <- TRUE
library(openxlsx)
devtools::load_all()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::READ_ALL
## Read sheets ####
# Also build specific !
tbls <- getSheetNames(file)
prefix <- c("v__","t__","i__","f__","m__","l__","e__","c__")
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::Stocks
if(verbose) cat("Read stocks\n") ### Stocks ####
namSp <- tbls[grepl("Stock__", tbls)]
nSp <- length(namSp)

# require tbls, file, namSp
Stock <- vector(mode = "list", length = nSp)
names(Stock) <- namSp
rm(nSp)

for(stock in namSp){
  if(verbose) cat(" ", stock)
  stock <- namSp[2] # TODO : dev to rm

  raw <- read.xlsx(file, sheet = stock, startRow = 35, rowNames=FALSE,
                   colNames=FALSE, skipEmptyRows = FALSE, skipEmptyCols = FALSE)
  raw[is.na(raw)] <- ""

  # TODO : only read recode for 1st sheet ? similare rec for all sheets ?
  rec <- read.xlsx(file, sheet = stock, rows = 4:38, cols = 1:4,
                   rowNames=FALSE, skipEmptyRows = FALSE, skipEmptyCols = FALSE)
  rec[is.na(rec)] <- ""
  rec[rec[,2] == "", 2] <- rec[rec[,2] == "", 1]
  vars <- rec$Variable # keeping name for later
  rec <- rec[rec$Variable != rec$Alias,] # rm identical name
  rec$Multi <- as.numeric(rec$Multi)

  # split tables with prefix
  indicRow <- apply(raw,1,function(x) any(substring(x,1,3) %in% prefix))
  indicTbl <- cumsum(apply(raw,1,function(x) all(x=="")))
  sepTabl <- unname(split(raw[indicRow,],indicTbl[indicRow]))
  # rm(raw, indicRow, indicTbl) # TODO : dev rm

  # format tables as correct list
  sepTabl <- lapply(sepTabl, function(x){
    x <- x[apply(x, 2, function(x) any(x != ""))]
    browser()
    if(substring(as.character(x[1,1]),1,3) %in% prefix){ # 1D
      x <- x[,1:((1:ncol(x))[!substring(x[1,],1,3) %in% prefix][1])]
      x <- IAM:::twoDto1D(x, "1D")
      x$value <- as.numeric(x$value)
    } else { # 2D
      x <- x[,apply(x,2,function(y) any(substring(y,1,3) %in% prefix))]
      x <- IAM:::twoDto1D(x, "2D")
    }
    browser()
    if(length(unique(x$v)) > 1){
      x <- lapply(split(x, x$v), function(d) { d$v <- NULL ; return(d) })
    } else {
      nme <- x$v[1]
      x <- list(x[,colnames(x) != "v", drop = FALSE])
      names(x) <- nme
    }

    return(x)
  })
  sepTabl <- unlist(sepTabl, recursive = FALSE)

  # Recode names
  nms <- substring(names(sepTabl), 4)
  nms[nms %in% rec$Alias] <- rec$Variable[match(nms[nms %in% rec$Alias],
                                                rec$Alias)]
  names(sepTabl) <- nms
  rm(nms) # TODO : dev rm
  # apply rec multiple
  for(n in names(sepTabl)){
    mult <- rec$Multi[rec$Variable == n]
    if (n %in% rec$Variable && mult > 1){
      sepTabl[[n]]$value <- sepTabl[[n]]$value * mult
    }
    tmp <- NA # TODO : dev rm
  }
  rm(n, tmp, mult, vars, rec) # TODO : dev rm


  Stock[[stock]] <- sepTabl
  rm(sepTabl, stock)
  if(verbose) cat(".")
}
if(verbose) cat("\n")
# return Stock
# tbls <- tbls[!grepl("^Stock__", tbls)]

# microbenchmark::microbenchmark(
# a = paste("FLWfmi_",as.vector(t(outer(paste("S",1:4,sep=""),paste("M",1:4,sep=""),paste,sep=""))),sep=""),
# ap = paste0("FLWfmi_",as.vector(t(outer(paste0("S",1:4),paste0("M",1:4),paste0)))),
# b = paste0("FLWfmi_",sort(outer(paste0("S",1:4),paste0("M",1:4),paste0))),
# times = 100
# )








# file <- "dev/data/inputFile.xlsx"
# t_init <- 2020
# t_hist_max <- 2020
# nbStep <-  5
# desc <- "Med tryhard"
# folderFleet = "dev/data/fleets"
# verbose <- TRUE
# out <- IAM:::read.input("dev/data/inputFile.xlsx",
#                         t_init = t_init, nbStep = nbStep,
#                         t_hist_max = t_hist_max, desc = desc,
#                         folderFleet = folderFleet, verbose = verbose)

