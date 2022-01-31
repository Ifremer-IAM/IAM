

devtools::load_all()
library(tidyverse)
library(viridis)
library(ggthemes)

setGestion_tac <- function(argum, species, delay = 1, bound = c(0,0), tac){

  argum@arguments$Gestion$active <- 1
  argum@arguments$Gestion$espece <- species
  argum@arguments$Gestion$delay <- delay
  argum@arguments$Gestion$sup <- bound[1]
  argum@arguments$Gestion$inf <- bound[2]
  argum@arguments$Gestion$tac <- tac

  return(argum)
}

load("dev/data/inputIFR.RData")
load("dev/data/argumIFR.RData")
argum1984 <- setGestion_tac(argum1984, "COR", 2, c(100, -100),
                            c(NA, 3000, 3000, 3001, 3000, 3000, 3000, 3000, 3000, 3000, 3000, 3000))

argum1984@arguments$Recruitment$ARC$typeMODsr <- "Hockey-Stick"
sim1984 <- IAM::IAM.model(objArgs = argum1984, objInput = input1984,
                          verbose = TRUE, force_t = 4)

# TACF ####

nbV <- length(input1984@specific$Fleet)
nbT <- input1984@specific$NbSteps

input1984@input$COR$Lref_f_e %>%
  as.data.frame() %>% mutate(perc = ./sum(.), Lref = sum(.)*perc) %>%
  select(-1) -> TACtmp

scenar <- c(1, 0.95, 0.93, 0.90, 1, 1, 1, 1, 1, 1, 1, 1, 1)
TACtmp <- matrix(TACtmp$Lref, ncol = nbT+1, nrow = nbV, byrow = FALSE)
for(y in 2:(nbT+1)){
  TACtmp[, y] <- TACtmp[, y-1] * scenar[y]
}
TACf <- TACtmp[, -1];  TACtot = colSums(TACtmp[, -1])
dimnames(TACf) <- list(input1984@specific$Fleet,input1984@specific$times)
TACf
rm(TACtmp, y, nbV, nbT, scenar)

devtools::load_all()
load("dev/data/inputIFR.RData")
load("dev/data/argumIFR.RData")
argum1984@arguments$Gestion$espece <- "COR"
argum1984@arguments$Gestion$inf <- -100
argum1984@arguments$Gestion$sup <- 100
argum1984@arguments$Gestion$active <- 1

TACf <- list(COR = TACf) ; TACtot <- list(COR = TACtot)

sim1984 <- IAM::IAM.model(objArgs = argum1984, objInput = input1984, TACbyF=TACf, TACtot = TACtot, updateE = 1,
                          verbose = TRUE, force_t = 4)

sim1984@output$effort1_f_m
sim1984@output$effort2_f_m

resF <- format_var("F", sim1984)  %>% filter(species == "COR")
  # filter(metier == "Filet_DP", fleet == "Antea", year >= 9) %>%
  # print(n = 22) %>%

resF %>%  #filter(metier == "Filet_COR", fleet == "Thalassa", year < 5) %>%
  filter( year < 5) %>%
ggplot(., aes(x=year, y=value, color = age)) +
  geom_line() +
  # facet_grid(fleet ~ metier ) +
  ggtitle("SSB with mean recrutment drawn in rnorm (sd = 1e6)")



dim(sim1984@outputSp$GVLcom_f_m_e$COR)

sim1984@output$gpmargin_f

sim1984@output$cs_f
sim1984@output$ps_f
sim1984@output$psAct_f
# Grouping variable for plot ####


# cas vide PQuot
format_varsp("PQuot", sim1984)
# cas simple SSB
format_varsp("SSB", sim1984)
# cas complexe F
format_varsp("F", sim1984) # %>% filter(is.na(value)) %>% head # TODO why NA here ?
format_varsp("N_S1M1", sim1984)

# demonstration de la puissance !!!!
format_vareco(name = "ratio_gp_K_f", sim1984)
format_vareco(name = "reconcilSPP", sim1984)

# group output and outputsp in one single large element of doom
# details : remove reconcilSPP because value is character and it's crap
simpl <- function(sim){
  output <- lapply(names(sim@output), format_vareco, sim) %>%
    do.call(rbind, .) %>%
    filter( variable != "reconcilSPP") %>%
    mutate( value = as.numeric(value),year = as.numeric(year) )

  outputsp <- lapply(names(sim@outputSp), format_varsp, sim) %>%
    do.call(rbind, .)%>%
    mutate( value = as.numeric(value),year = as.numeric(year) )

  bind_rows(output, outputsp)
}


# Test when variable does not exist return NULL
format_var("MLA",sim1984)

# ici on recup toutes les variables pour une simul
# permet d'imaginer un rbind de plusieurs simuls
res <- lapply(names(sim1984@outputSp), format_var, sim1984) %>%
  do.call(rbind, .)

# ici on filtre une simul pour plot selon metier flotte.
res %>% filter(variable == "GVL_f_m_e") %>%
  mutate( year = as.numeric(year) ) %>%
  mutate( value = as.numeric(value) ) %>%
  ggplot(., aes(x=year, y=value, color = fleet)) +
  geom_line() +
  facet_grid(species ~ metier)

# et la on peut meme trier facilement,
# enlever les colonnes inutiles et mettre au bon format
# le format de base etant celui pour ggplot
format_var(name = "allocEff_f_m", sim1984) %>%
  filter( fleet == "Alis") %>%
  discard( ~n_distinct(.) == 1) %>%
  pivot_wider(names_from = year, values_from = value)





# Testing replication ####
devtools::load_all()
load("dev/data/inputIFR.RData")
load("dev/data/argumIFR.RData")

NITER <- 10
XSAsp <- names(argum1984@arguments$Recruitment)

pb <- txtProgressBar(min = 0, max = NITER, style = 3)
simuls <- lapply(1:NITER,function(k) {
  argum <- argum1984
  for(nms in XSAsp){
    val <- argum1984@arguments$Recruitment[[nms]]$parAmodSR
    argum@arguments$Recruitment[[nms]]$parAmodSR <- rnorm(1, val, 1e6)
  }
  setTxtProgressBar(pb, k)
  return(IAM.model(objArgs = argum, objInput = input1984))
})
close(pb)

# take 2 minutes...not a good idea.
pb <- txtProgressBar(min = 0, max = NITER, style = 3)
res <- lapply(1:NITER,function(k) {
  setTxtProgressBar(pb, k)
  format_var("SSB", simuls[[k]]) %>%
    mutate( value = as.numeric(value),year = as.numeric(year) ) %>%
    add_column(sim = k, .before = "variable")
  }) %>%
  do.call(rbind, .) %>%
  filter(species != "DAR")
close(pb)

# https://stackoverflow.com/questions/34749859/ggplot2-shading-envelope-of-time-series
condquant <- res %>%
  discard( ~n_distinct(.) == 1) %>%
  group_by(year, species) %>%
  do(quant = quantile(.$value, probs = seq(0,1,.25)), probs = seq(0,1,.25)) %>%
  unnest(cols = c(quant, probs)) %>%
  mutate(delta = 2*round(abs(.5-probs)*100)) %>%
  group_by(year, delta, species) %>%
  summarize(quantmin = min(quant), quantmax= max(quant), quantmean  = mean(quant)) %>%
  add_column(value = 1)

res  %>% filter(sim == 500) %>%
  ggplot(., aes(x=year, y=value)) +
  geom_ribbon(data = condquant, aes(x = year, ymin = quantmin, ymax = quantmax,
                                    group = reorder(delta, -delta), fill = as.numeric(delta)),
              alpha = .5) +
  scale_fill_continuous( limits = c(50, 100), name = "Density") +
  geom_line() +
  facet_grid(species ~ .) + ggtitle("SSB with mean recrutment drawn in rnorm (sd = 1e6)")
  #+ theme_dark() +
  # theme(panel.background = element_rect(fill = "steelblue"))


