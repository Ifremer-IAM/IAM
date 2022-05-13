#' Scrip for testing the IAM.test_plot function
#'
#' author : maxime Jaunatre 29/03/2022
#'
#' script of the IAM package.
#'
devtools::load_all() # Load IAM
data("IAM_input_2009")
data("IAM_argum_2009")
sim2009 <- IAM::IAM.model(objArgs = IAM_argum_2009, objInput = IAM_input_2009,
                          verbose = FALSE, force_t = NULL)
sim2009 <- IAM.format(sim2009, "summary",sim_name = "scenar1A", n= 2)
sim2009 <- IAM.format_quant(sim2009) %>%
  dplyr::filter(variable != "effort2_f")

options(shiny.launch.browser = .rs.invokeShinyWindowExternal)
plot <- IAM.test_plot(sim2009)
