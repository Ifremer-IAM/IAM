#' @import shiny
#' @rawNamespace import(shinyjs, except = c(runExample, alert)) #https://github.com/daattali/shinyjs/issues/127
#' @importFrom shinyWidgets pickerInput useSweetAlert
#'
#' @noRd
mod_BioplotUI <- function(id){
  ns <- NS(id)
  tagList(
    useSweetAlert(),
    hidden(
      pickerInput(
        inputId = ns("variable"), label = "Variables",
        choices = "Number", selected = "Number", multiple = TRUE,
        options = list(`actions-box` = TRUE)
      ),
      pickerInput(
        inputId = ns("species"), label = "Species",
        choices = "Fish", selected = "Fish"
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
#' Permet de faire les plots de Fbar, SSB et L_et en fonction du temps.
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
mod_Bioplot_serv <- function(id, x){

  moduleServer(
    id,
    ## Below is the module function
    function(input, output, session) {

      def <- reactiveValues(sim_name = NULL, species = NULL,
                            fullvar = NULL, variable = NULL)

      observe({
        def$sim_name <- x$sim_name
      })

      observeEvent(x$bioeco, {
        toggle("species", condition = x$bioeco == "Biologic")
        toggle("variable", condition = x$bioeco == "Biologic")
      })

      observe({
        var <- get_golem_options("input")

        for (axe in c("variable", "sim_name", "species")) {
          if (all(is.na(var[[axe]])) ) {
            tmp <- switch(axe,
                          species = "Fish",
                          sim_name = "Scenar",
                          variable = "Number")
          } else {
            tmp <- as.character(unique(var[[axe]]))
            tmp <- tmp[!is.na(tmp)]

            if(axe == "variable"){
              tmp = get_golem_options("Biovars")
            }
            if(axe == "species" & !is.null(get_golem_options("dinsp"))){
              tmp <- get_golem_options("dinsp")
            }
          }

          if(axe == "species"){n_sel = 1} else {n_sel = seq_along(tmp)}
          updatePickerInput(
            session, inputId = axe, choices = tmp, selected = tmp[n_sel]
          )

          def[[axe]] <- tmp[n_sel]
          toggle(axe, condition = !all(is.na(var[[axe]])))
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

      observeEvent(input$variable, {def$variable <- input$variable })
      observeEvent(input$species, {def$species <- input$species })


      observe({
        N_vr <- !is.null(def$variable)
        N_sp <- !is.null(def$species)
        N_sn <- !is.null(def$sim_name)

        if (any(def$species == "Fish")) def$species <- NULL
        if (any(def$sim_name == "Apocalypse")) def$sim_name <- NULL
        if (any(def$variable == "Number")) def$variable <- NULL

        if (is.null(def$fullvar)) {
          df <- tibble()
        } else {
          df <- def$fullvar %>%
            dplyr::filter(if(N_sp) .data$species %in% def$species else TRUE) %>%
            dplyr::filter(if(N_sn) .data$sim_name %in% def$sim_name else TRUE) %>%
            dplyr::filter(if(N_vr) .data$variable %in% def$variable else TRUE) %>%
            dplyr::filter(.data$year >= x$time[1], .data$year <= x$time[2])
        }

        x$var <- df # for export
        # print(df)

        col = ifelse(x$colors, "white", "gray")

        if (nrow(df) > 0) {
          # browser()
          p <- ggplot(df, aes(x = .data$year, y = .data$median)) +
            { if (N_sn) facet_grid(.data$variable ~ .data$sim_name,
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

        x$bio_plot <- p
      })
    }
  )
}


