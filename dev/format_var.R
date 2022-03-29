
library(IAM)
devtools::load_all()

# testing the full process !
library(IAM)
a <- Sys.time()
IAM_argum_1984@arguments$Recruitment$COR$wnNOISEmodSR <- 5
nsim_statu_quo <- sapply(1:5, function(x){
  tmp <- IAM::IAM.model(objArgs = IAM_argum_1984, objInput = IAM_input_1984)
  tmp@outputSp$SSB$COR <- tmp@outputSp$SSB$COR + rnorm(12, sd = 500) # just to test the ribbon
  tmp <- IAM.format(tmp, "summary", "statuquo", x)
  return(tmp)
}, simplify = FALSE
)
nsim_tabl <- do.call(rbind, nsim_statu_quo)
nsim_tabl <- IAM.format_quant(nsim_tabl, probs = c(0.05, 0.95))
b <- Sys.time()
b - a # 1.17 min

nsim_tabl %>%
  filter(variable == "SSB", species != "DAR") %>%
ggplot(aes(x = year, y = median)) +
  geom_ribbon(aes(ymin = quant1, ymax = quant2), fill = "lightblue")+
  geom_line() +
  geom_line(aes(y = value), linetype = 2) +
  facet_grid(species ~ sim_name) +
  NULL



