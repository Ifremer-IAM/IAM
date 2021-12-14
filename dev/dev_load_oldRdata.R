
#' Script to load old .Rdata created with lost packages.
#' in "IAMXXX" you need to replace XXX with the version of the package.

rm(list = ls()) ; .rs.restartR() # clear everything
library(devtools)
tst <- namespace::makeNamespace("IAM65") # create an empty namespace to load S4
assign("H2G2", "DON'T PANIC", env = tst)   # require interactive run !
base::namespaceExport(tst, ls(tst))
IAM65::H2G2

# library(IAMtest64)
load("C:/Users/mjaunatr/Desktop/Docs passation/inputs/inputs_GG/input2016_newStocks_MNZstatic.RData")
# I can attribute S4 class to another package. This should be impossible to save my sanity.
attributes(class(input2016))$package <- "IAM"
