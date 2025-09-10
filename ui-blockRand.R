# ui-blockRand.R

fluidPage(
  includeCSS("www/style/theme.css"),
  fluidRow(
    column(
      width = 3,
      box(
        title = tagList(icon("cloud-upload"), "Sample Information"),
        status = "primary",
        solidHeader = TRUE,
        width = NULL,
        fileInput(
          "uploadData",
          "Upload Sample Information",
          accept = c(".csv")
        ),
        helpText("The file needs to be saved in .csv format. The first column should be the sample ID, with factors to be considered starting from the second column."),
        a(href = "data/example_input_for_block_randomization.csv",
          "Download example sample information",
          download = NA,
          target = "_blank"
        )
      ),
      box(
        title = tagList(icon("cogs"), "Experiment Settings"),
        status = "primary",
        width = NULL,
        solidHeader = TRUE,
        selectInput("experiment_type", "Experiment Type", choices = c("Label-free", "TMT-labeling")),
        numericInput(inputId = "batch_volume", label = "Number of all samples in a batch", value = 0),
        numericInput("n_IRs", "No. of IR (Internal Reference) samples in a batch", value = 0),
        sliderInput("n_iterations", "Optimization level", value = 100, min = 50, max = 1000),
        do.call(actionBttn, c(
          list(
            inputId = "run",
            label = "Run Experiment Design",
            icon = icon("play")
          ),
          actionBttnParams
        ))
      )
    ),
    column(
      width = 9,
      fluidRow(
        column(
          width = 8,
          box(
            title = tagList(icon("table"), "Sample Information Table"),
            width = NULL,
            solidHeader = TRUE,
            status = "primary",
            uiOutput("pre_inputData")
          )),
        column(
          width = 4,
          box(
            title = tagList(icon("bar-chart"), "Factor Distributions"),
            width = NULL,
            solidHeader = TRUE,
            status = "primary",
            uiOutput("factorDistribution")
          )
        )
      ),
      box(
        title = tagList(icon("table"), "Batch and Channel Assignation"),
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        uiOutput("pre_sample_batch"),
        downloadButton("download_all", label = "Download all results")
      ),
      box(
        title = tagList(icon("table"), "Batch Design Matrix"),
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        uiOutput("batch_design")
      )
    )
  )
)
