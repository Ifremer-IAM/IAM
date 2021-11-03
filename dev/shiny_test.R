library(shiny)
library(shiny.i18n)
library(shinyWidgets)
library(shinythemes)


a <- c(
  "A1" = "A1",
  "A2" = "A2",
  "A3" = "A3",
  "A4" = "A4"
)

sp <- c("HKE" = "HKE", "MUT" = "MUT")

# Define UI
ui <- fluidPage(
  fluidRow(
    column(
      width = 3,
      div(
        class = "option-group",
        radioButtons(
          "recrut", "Recruitement",
          choices = c("SP1", "SP2", "SP3"), inline = TRUE
        ),
        radioButtons(
          "noise", "Noise dist.",
          choices = c("Norm", "LogN"), inline = TRUE
        ),
        radioButtons(
          "plot_type", "Plot type",
          c("base", "ggplot2"),
          inline = TRUE
        ),
        conditionalPanel(
          "input.plot_type === 'base'",
          selectInput("plot_scaletype", "Scale type",
            c(
              "normal" = "normal",
              "log" = "log",
              "x factor" = "x_factor",
              "datetime" = "datetime"
            ),
            selectize = FALSE
          )
        ),
        conditionalPanel(
          "input.plot_type === 'ggplot2'",
          selectInput(
            "ggplot_scaletype", "Scale type",
            c(
              "normal" = "normal",
              "reverse (scale_*_reverse())" = "reverse",
              "log10 (scale_*_log10())" = "log10",
              "log2 (scale_*_continuous( trans=log2_trans()))" = "log2",
              "log10 (coord_trans())" = "log10_trans",
              "log2 (coord_trans())" = "log2_trans",
              "coord_cartesian()" = "coord_cartesian",
              "coord_flip()" = "coord_flip",
              "coord_fixed()" = "coord_fixed",
              "coord_polar() (doesn't work)" = "coord_polar",
              "x factor" = "x_factor",
              "date and time" = "datetime"
            ),
            selectize = FALSE
          ),
          selectInput(
            "ggplot_facet", "Facet",
            c(
              "none" = "none",
              "wrap" = "wrap",
              "grid x" = "grid_x",
              "grid y" = "grid_y",
              "grid xy" = "grid_xy",
              "grid xy free" = "grid_xy_free"
            ),
            selectize = FALSE
          )
        )
      ),
      div(
        class = "option-group",
        checkboxInput(
          "iter", "Iterative",
          value = FALSE
        ),
        conditionalPanel(
          "input.iter",
          selectInput("n_iter", "Number of iterations",
            c(
              "dix" = "10",
              "cent" = "100",
              "mille" = "1000"
            ),
            selectize = FALSE
          ),
          selectInput("variables", "Output variables",
            c(
              "normal" = "B",
              "r" = "SSB",
              "log10" = "log10",
              "log2" = "Ctot",
              "log10)" = "ytot"
            ),
            selectize = FALSE
          )
        )
      )
    ),
    column(
      width = 6,
      fluidRow(
        column(
          width = 3,
          div(
            class = "option-group",
            checkboxInput("manag", "Management", value = FALSE),
            conditionalPanel(
              "input.iter",
              selectInput("n_iter", "Number of iterations",
                c(
                  "dix" = "10",
                  "cent" = "100",
                  "mille" = "1000"
                ),
                selectize = FALSE
              ),
              selectInput("variables", "Output variables",
                c(
                  "normal" = "B",
                  "r" = "SSB",
                  "log10" = "log10",
                  "log2" = "Ctot",
                  "log10)" = "ytot"
                ),
                selectize = FALSE
              )
            )
          )
        ),
        column(
          width = 3,
          div(
            class = "option-group",
            checkboxInput("eco", "Economic", value = FALSE),

            # conditionalPanel("input.eco",
            radioGroupButtons(
              inputId = "Id066",
              label = "Label",
              choices = c(
                "Complete",
                "DCF"
              ),
              justified = TRUE
            ),
            awesomeRadio(
              inputId = "cc",
              label = "Adj",
              choices = c(
                "1",
                "2"
              ),
              selected = "1",
              inline = TRUE
            ),
            awesomeRadio(
              inputId = "ss",
              label = "ue_choice",
              choices = c(
                "1",
                "2"
              ),
              selected = "1",
              inline = TRUE
            ),
            awesomeRadio(
              inputId = "ii",
              label = "oths",
              choices = c("0", "1"),
              selected = "1",
              inline = TRUE
            ),
            awesomeRadio(
              inputId = "kk",
              label = "othsFM",
              choices = c("0", "1"),
              selected = "1",
              inline = TRUE
            )
            # )
          )
        )
      ),
      div(
        class = "option-group",
        checkboxInput("scenar", "Scenario", value = FALSE),

        # conditionalPanel("input.scenar",
        selectInput("variables", "Output variables",
          a,
          selectize = FALSE
        )
        # )
      )
    ),
  )
)

# Define server function
server <- function(input, output) {

}

# Create Shiny object
shinyApp(ui = ui, server = server)
