

devtools::load_all()
library(tidyverse)
library(viridis)
library(ggthemes)

# TACF ####
devtools::load_all()
data("IAM_argum_2009")
data("IAM_input_2009")

nbV <- length(IAM_input_2009@specific$Fleet)
nbT <- IAM_input_2009@specific$NbSteps

Lref <- IAM_input_2009@input$COR$Lref_f_e
percLref <- Lref / sum(Lref)

scenar <- c(1, 0.95, 0.93, 0.90, 1, 1, 1, 1, 1, 1, 1, 1, 1)
TACtmp <- matrix(Lref, ncol = nbT+1, nrow = nbV, byrow = FALSE)
for(y in 2:(nbT+1)){
  TACtmp[, y] <- TACtmp[, y-1] * scenar[y]
}
TACf <- TACtmp[, -1];  TACtot = colSums(TACtmp[, -1])
dimnames(TACf) <- list(IAM_input_2009@specific$Fleet,IAM_input_2009@specific$times)
TACf
rm(TACtmp, y, nbV, nbT, scenar)


IAM_argum_2009@arguments$Gestion$espece <- "COR"
IAM_argum_2009@arguments$Gestion$inf <- -100
IAM_argum_2009@arguments$Gestion$sup <- 100
IAM_argum_2009@arguments$Gestion$active <- 1

TACf <- list(COR = TACf) ; TACtot <- list(COR = TACtot)

sim2009 <- IAM::IAM.model(objArgs = IAM_argum_2009, objInput = IAM_input_2009,
                          TACbyF=TACf, TACtot = TACtot, updateE = 1,
                          verbose = TRUE)

sim2009@output$effort1_f_m
sim2009@output$effort2_f_m

resF <- IAM.format(sim2009, "F")  %>% filter(species == "COR")
  # filter(metier == "Filet_DP", fleet == "Antea", year >= 9) %>%
  # print(n = 22) %>%

resF %>%  #filter(metier == "Filet_COR", fleet == "Thalassa", year < 5) %>%
  filter( year < 5) %>%
ggplot(., aes(x=year, y=value, color = age)) +
  geom_line() +
  # facet_grid(fleet ~ metier ) +
  ggtitle("SSB with mean recrutment drawn in rnorm (sd = 1e6)")


dim(sim2009@outputSp$GVLcom_f_m_e$COR)

sim2009@output$gpmargin_f

sim2009@output$cs_f
sim2009@output$ps_f
sim2009@output$psAct_f

