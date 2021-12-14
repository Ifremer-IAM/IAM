file <- "dev/data/inputFile.xlsx"
t_init <- 2020
t_hist_max <- 2020
nbStep <-  5
desc <- "Med tryhard"
folderFleet = "dev/data/fleets"
verbose <- TRUE
library(openxlsx)

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
for(stock in namSp){

  stock <- namSp[1] # TODO : dev to rm
  raw <- read.xlsx(file, sheet = stock, startRow = 35, rowNames=FALSE, colNames=FALSE,
                              skipEmptyRows = FALSE, skipEmptyCols = FALSE)
  raw[is.na(raw)] <- ""
  # TODO : only read recode for 1st sheet ? similare recode for all sheets ?
  recode <- read.xlsx(file, sheet = stock, rows = 5:38, cols = 1:4,
                      rowNames=FALSE, skipEmptyRows = FALSE, skipEmptyCols = FALSE)
  recode[is.na(recode)] <- ""
  recode[recode[,2] == "", 2] <- recode[recode[,2] == "", 1]
  recode <- recode[recode$Multi != "1",] # we don't need to apply multiple of 1...
  recode$Multi <- as.numeric(recode$Multi)

  # split tables with prefix
  indicRow <- apply(raw,1,function(x) any(substring(x,1,3) %in% prefix))
  indicTbl <- cumsum(apply(raw,1,function(x) all(x=="")))
  sepTabl <- unname(split(raw[indicRow,],indicTbl[indicRow]))

  sepTabl <- lapply(sepTabl, function(x){
    x <- x[apply(x, 2, function(x) any(x != ""))]

    if(substring(as.character(x[1,1]),1,3) %in% prefix){ # 1D
      x <- x[,1:((1:ncol(x))[!substring(x[1,],1,3) %in% prefix][1])]
      x <- IAM:::twoDto1D(x, "1D")
    } else { # 2D
      x <- x[,apply(x,2,function(y) any(substring(y,1,3) %in% prefix))]
      x <- IAM:::twoDto1D(x, "2D")
    }
    return(x)
  })




  Stock[[stock]] <- raw
}
# return Stock
tbls <- tbls[!grepl("^Stock__", tbls)]

microbenchmark::microbenchmark(
a = paste("FLWfmi_",as.vector(t(outer(paste("S",1:4,sep=""),paste("M",1:4,sep=""),paste,sep=""))),sep=""),
ap = paste0("FLWfmi_",as.vector(t(outer(paste0("S",1:4),paste0("M",1:4),paste0)))),
b = paste0("FLWfmi_",sort(outer(paste0("S",1:4),paste0("M",1:4),paste0))),
times = 100
)
