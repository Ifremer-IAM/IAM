#'
#' Script to load old .Rdata created with lost packages.
#' in "IAMXXX" you need to replace XXX with the version of the package.
#'
#' Maxime Jaunatre
#'
.rs.restartR() # clear everything
#'
#' You need to replace "IAM64 with the name of the missing package"
#'
tst <- namespace::makeNamespace("IAM64") # create an empty namespace to load S4
base::assign("H2G2", "DON'T PANIC", env = tst)   # require interactive run !
base::namespaceExport(tst, ls(tst))
IAM64::H2G2

load(paste0("C:/Users/mjaunatr/Desktop/Docs passation/inputs/",
            "inputs_GG/input2016_newStocks_MNZstatic.RData"))
# I can attribute S4 class to another package.
# This should be impossible to save my sanity.
attributes(class(input2016))$package <- "IAM"
attributes(class(argum2016))$package <- "IAM"
