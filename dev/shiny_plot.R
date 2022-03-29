#' Scrip for testing the IAM.test_plot function
#'
#' author : maxime Jaunatre 29/03/2022
#'
#' script of the IAM package.


devtools::load_all() # Load IAM
data("IAM_input_1984")
data("IAM_argum_1984")
sim1984 <- IAM::IAM.model(objArgs = IAM_argum_1984, objInput = IAM_input_1984,
                          verbose = FALSE, force_t = NULL)
sim1984 <- IAM.format(sim1984, "summary",sim_name = "scenar1A", n= 2)
sim1984 <- IAM.format_quant(sim1984) %>%
  dplyr::filter(variable != "effort2_f")

options(shiny.launch.browser = .rs.invokeShinyWindowExternal)
plot <- IAM.test_plot(sim1984)
