fluidPage(
  includeCSS("www/style/theme.css"),
  fluidRow(
    column(
      3,
      box(
        title = tagList(icon("cogs"), "Differential Expression Methods"),
        id = "DEselect",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        selectInput(
          "diffMethod",
          "Differential Expression Method",
          choices = c("LIMMA" = "limma")
        ),
        uiOutput("groups2"),
        uiOutput("compSelect"),
        uiOutput("group3"),
        uiOutput("selectInputContainer"),
        radioButtons(
          inputId = "ImputationSelection",
          label = "Missing Value Imputation",
          choices = c(
            "No Imputation" = "NA",
            "Imputation" = "Imputation",
            "Data without NAs" = "CleanData"
          ),
          selected = "NA"
        ),
        conditionalPanel(
          condition = "input.ImputationSelection == 'Imputation'",
          numericInput(
            inputId = "filterValue",
            label = "Enter Filter Value",
            value = 1
          )
        ),
        # Popovers for each radio button choice
        radioTooltip("ImputationSelection", choice = "NA", title = "No imputation on whole data"),
        radioTooltip("ImputationSelection", choice = "Imputation", title = "See documentation for details"),
        radioTooltip("ImputationSelection", choice = "CleanData", title = "No imputation and perform DE only on data without NAs"),
        do.call(
          actionBttn, c(
            list(
              inputId = "DETestType",
              label = "Run Differential Expression Test",
              icon = icon("play")
            ),
            actionBttnParams
          )
        )
      )
    ),
    column(
      9,
      box(
        title = tagList(icon("table"), "Differential Expression Result Table"),
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        uiOutput("DEResultTable")
      ),
      tabBox(
        title = "",
        width = NULL,
        tabPanel(
          title = tagList(icon("bar-chart"), "Moving SD"),
          uiOutput("DE_distributionUI")
        ),
        tabPanel(
          title = tagList(icon("bar-chart"), "SD within group"),
          uiOutput("DE_OldSDUI")
        ),
        tabPanel(
          title = tagList(icon("bar-chart"), "Volcano Plot"),
          uiOutput("DE_volcanoPlot")
        ),
        tabPanel(
          title = tagList(icon("sitemap"), "Heatmap"),
          uiOutput("DE_dendUI")
        )
      )
    )
  )
)
