
# Copy from https://ropensci.org/blog/2019/12/08/precompute-vignettes/
library(knitr)
unlink(c("vignettes/figure"), recursive = TRUE) # Remove the directories if they exist
dir.create("vignettes/figure")
knit("vignettes/Target_Fmsy.Rmd.orig", "vignettes/Target_Fmsy.Rmd")
file.copy(list.files("figure", full.names = TRUE), to = "vignettes/figure")
unlink(c("figure"), recursive = TRUE)
