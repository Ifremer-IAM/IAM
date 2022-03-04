library(shiny)
library(shinyWidgets)
library(rhandsontable)
library(shinyjs)


# Plotting Rec and/or hide useless parameters ####
# 1 -> recrutement constant moyen (rec~a)
# 2 -> Hockey stick (rec ~ (si (ssb<=b) a*ssb sinon a*b))
# 3 -> Beverton & Holt (rec ~ a*ssb/(b+ssb))
# 4 -> Ricker (rec ~ a*ssb*exp(-b*ssb))
# 5 -> Shepherd (rec ~ a*ssb/(1+ (ssb/b)^c))
# 6 -> Hockey Stick Quadratic (rec ~ (si (ssb<=b*(1-c)) a*ssb ; si (b*(1-c)<ssb<b*(1+c)) a*(ssb-((ssb-b*(1-c))^2)/(4*b*c)) ; sinon a*b))
# 7 -> Hockey Stick Smooth (rec ~ a*(ssb+sqrt(b^2+g)-sqrt((ssb-b)^2+g)), avec g=0.001 )


# just to test it for now
AllVarRep = c(
  "B", "SSB", "Ctot", "Ytot", "Yfmi", "Ffmi", "Zeit", "Fbar", "Foth",
  "mu_nbds", "mu_nbv", "N", "Ystat", "Lstat", "Dstat", "Eff",
  "GVL_fme", "StatGVL_fme", "GVLtot_fm", "GVLav_f", "vcst_fm", "vcst_f",
  "rtbs_f", "gp_f", "ps_f", "gcf_f", "gva_f", "cs_f", "sts_f", "rtbsAct_f",
  "csAct_f", "gvaAct_f", "gcfAct_f", "psAct_f", "stsAct_f", "ccwCr_f",
  "GVLtot_f", "wagen_f", "L_efmit", "D_efmit",
  "Fr_fmi", "C_efmit", "P", "Pstat"
)

library(IAM)
data("IAM_argum_1984")
# cas avec les choix de tout et n'importe quoi
IAM_argum_1984@arguments$Scenario$ALLscenario <- c("Scenario1", "Scenario2")
IAM_argum_1984 <- IAM.editArgs_Scenar(IAM_argum_1984, 2)

IAM_argum_1984 <- IAM.editArgs_Gest(
  IAM_argum_1984, active = TRUE, control = "Nb trips", target = "Fbar",
  espece = "DAR", delay = 6, type = "x", update = TRUE, bounds = c(100, 1))

IAM_argum_1984 <- IAM.editArgs_Eco(IAM_argum_1984, dr = 0.04, perscCalc = 3)

IAM_argum_1984@arguments$Replicates$nbIter <- 300
IAM_argum_1984@arguments$Replicates$active <- 1

IAM_arg_app(IAM_argum_1984, AllVarRep)
summary(IAM_argum_1984)

IAM::IAM.args(IAM_argum_1984)

i <- IAM.args(IAM_argum_1984)
IAM.args(IAM_input_1984)


# enveloppe de loi normale.
library(tidyverse)
sd = 1
d <- data.frame(a = seq(10,11,length.out = 10)  ) %>%
  mutate(d = qnorm(0.05,mean=a,sd=sd)) %>%
  mutate(b = qnorm(0.95,mean=a,sd=sd)) %>%
  mutate(dl = qlnorm(0.05,mean=a,sd=sd)) %>%
  mutate(bl = qlnorm(0.95,mean=a,sd=sd))


d %>%
  ggplot(aes(x = a, y = a)) +
  # geom_ribbon(aes(ymin = b, ymax = d), fill = "azure3") +
  geom_ribbon(aes(ymin = a, ymax = bl), fill = "azure3") +
  geom_line() +
  NULL


