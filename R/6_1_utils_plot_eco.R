#' @import shiny
#' @rawNamespace import(shinyjs, except = c(runExample, alert)) #https://github.com/daattali/shinyjs/issues/127
#' @importFrom shinyWidgets pickerInput useSweetAlert
#'
#' @noRd
mod_EcoplotUI <- function(id){
  ns <- NS(id)
  tagList(
    hidden(
      pickerInput(
        inputId = ns("variable"), label = "Variables",
        choices = "Number"
      ),
      pickerInput(
        inputId = ns("fleet"), label = "Fleet",
        choices = "Boat", selected = "Boat", multiple = TRUE,
        options = list(`actions-box` = TRUE)
      ),
      actionBttn(
        inputId = "missing", label = "Missing variable",
        style = "pill", color = "warning"
      )
    )
  )
}

#' Module shiny pour la partie Bio de l'explo graphique.
#'
#' Permet de faire le plot au choix de
#'  nbv
#'  effort2
#'  GVLav
#'  gva
#'  gp
#'  wageg
#'  wagen
#' si une variable est entieremenet compose de NA, elle ne sera pas dispo
#'
#' Depend du petit x reactif qui contient:
#'  bioeco : tri sur le fait d'utiliser le module bio ou eco.
#'  sim_name : tri sur les scenarios a presenter
#'  checkribbon : option de representation des effets stockastiques
#'  time : tri sur la periode de simulation a representer
#'
#'  il contient aussi des variables de retour comme
#'  var : permet d'exporter le tableau trie
#'  bio_plot : permet de representer le plot
#'  eco_plot : pas affecte par cette fonciton
#'
#' @author Maxime Jaunatre
#' @import ggplot2
#'
#' @noRd
mod_Ecoplot_serv <- function(id, x){

  moduleServer(
    id,
    ## Below is the module function
    function(input, output, session) {

      def <- reactiveValues(sim_name = NULL, fleet = NULL,
                            fullvar = NULL, variable = NULL)

      observe({
        def$sim_name <- x$sim_name
      })

      observeEvent(x$bioeco, {
        toggle("fleet", condition = x$bioeco == "Economic")
        toggle("variable", condition = x$bioeco == "Economic")
      })

      observe({
        var <- get_golem_options("input")

        for (axe in c("variable", "sim_name", "fleet")) {
          if (all(is.na(var[[axe]])) ) {
            tmp <- switch(axe,
                          fleet = "Boat",
                          sim_name = "Scenar",
                          variable = "Number")
          } else {
            tmp <- as.character(unique(var[[axe]]))
            tmp <- tmp[!is.na(tmp)]
            if(axe == "variable"){
              tmp = get_golem_options("Ecovars")

              tmp <- var %>% dplyr::filter(.data$variable %in% tmp) %>%
                group_by(.data$variable) %>%
                summarise(all_na = all(is.na(.data$value))) %>%
                dplyr::filter(! .data$all_na) %>% pull(.data$variable)
            }
          }

          if(axe == "variable"){ n_sel = 1 } else { n_sel = seq_along(tmp) }
          updatePickerInput(
            session, inputId = axe, choices = tmp, selected = tmp[n_sel]
          )

          def[[axe]] <- tmp[n_sel]
          toggle(axe, condition = !all(is.na(var[[axe]])) &
                   x$bioeco == "Economic")
        }

        toggle("missing", condition = nrow(var) == 0)
        def$fullvar <- var

      })

      observeEvent(input$missing, {
        sendSweetAlert(
          session = session,
          title = "Missing variable",
          text = paste(
            "No variable was not computed by this simulation or",
            "lack every dimension used by this function."
          ),
          type = "Warning"
        )
      })

      observeEvent(input$variable, {def$variable <- input$variable})
      observeEvent(input$fleet, {def$fleet <- input$fleet })


      observe({
        N_vr <- !is.null(def$variable)
        N_fl <- !is.null(def$fleet)
        N_sn <- !is.null(def$sim_name)

        if (any(def$fleet == "Boat")) def$fleet <- NULL
        if (any(def$sim_name == "Apocalypse")) def$sim_name <- NULL
        if (any(def$variable == "Number")) def$variable <- NULL

        if (is.null(def$fullvar)) {
          df <- tibble()
        } else {
          df <- def$fullvar %>%
            dplyr::filter(if (N_fl) .data$fleet %in% def$fleet else TRUE) %>%
            dplyr::filter(if (N_sn) .data$sim_name %in% def$sim_name else TRUE) %>%
            dplyr::filter(if (N_vr) .data$variable %in% def$variable else TRUE) %>%
            dplyr::filter(.data$year >= x$time[1], .data$year <= x$time[2])
        }

        x$var <- df # for export
        # print(df)

        col = ifelse(x$colors, "white", "gray")

        if (nrow(df) > 0) {

          # Making the plot in question
          p <- ggplot(df, aes(x = .data$year, y = .data$median)) +
            { if (N_sn) facet_grid(.data$fleet ~ .data$sim_name,
                                   scales = "free_y") } +
            { if (x$chkribbon) geom_ribbon(aes(ymin = .data$quant1,
                                               ymax = .data$quant2),
                                           fill = col, alpha = .4) } +
            geom_line(size = .5) +
            geom_point(fill=col, size = .5) +
            geom_line(aes(y = .data$value), linetype = "dashed") +
            guides(x = guide_axis(angle = 90)) +
            IAM_theme(blue = x$colors) +
            NULL
        } else {
          p <- ggplot(NULL)
        }

        x$eco_plot <- p
      })
    }
  )
}


