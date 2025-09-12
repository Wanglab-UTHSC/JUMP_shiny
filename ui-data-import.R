# ui-data-import.R

fluidPage(
  includeCSS("www/style/theme.css"),
  fluidRow(column(
  3,
  #tabBox(
  box(
      title = tagList(icon("cloud-upload"), "Upload"),
      solidHeader = TRUE,
      status = "primary",
      width = NULL,
      radioButtons(
        "dataType",
        "Choose your file type",
        choices = c("jumpshiny",
                    "jumpq",
                    "jump_batch")
      ),
      radioButtons(
        "intensityType",
        "Type of intensity values",
        choices = c("raw", "log2"), 
        inline = T
      ),
      fileInput(
        "uploadExpressionData",
        "Upload Data",
        accept = c("text/csv",
                   "text/comma-separated-values,text/plain",
                   ".csv"),
        buttonLabel = "Upload...",
        placeholder = "No file has been uploaded."
      ),
      helpText("Text file in .tsv/.csv format, and the first column should be Accession Number, second column should be Gene Names, Third should be Description."),
      a(href = "data/example_data.zip",
        "Download example input expression table",
        download = NA,
        target = "_blank"
      )

  ),
  box(
    title = tagList(icon("tags"), "Meta information"),
    solidHeader = TRUE,
    status = "primary",
    width = NULL,
    fileInput(
      "uploadGroup",
      "Upload sample information",
      accept = c("text/csv",
                 "text/comma-separated-values,text/plain",
                 ".csv"),
      buttonLabel = "Upload...",
      placeholder = "No file has been uploaded."
    ),
    helpText("Text file in .tsv/.csv format, and the first column should be named as Sample."),
    # a(href = "data/meta_DIA.csv",
      # "Download example sample information table",
      # download = NA,
      # target = "_blank"
    # ),
    tags$br(),
    do.call(actionBttn, c(
      list(
        inputId = "confirmedGroupList",
        label = "Assign group information",
        icon = icon("play")),
      actionBttnParams
      )
    ),
    footer = helpText("JUMP Shiny expect first column to be sample name, second column to be group information, and so on.")
  ),
  box(
    title = tagList(icon("layer-group"), "Group selection"),
    solidHeader = TRUE,
    status = "primary",
    width = NULL,
    uiOutput("groupSelection")
  ),
  box(
    title = tagList(icon("info-circle"), "Data info"),
    solidHeader = TRUE,
    status = "info",
    width = NULL,
    uiOutput("importDataSummary")
  )
),
column(
  9,
  box(
    title = tagList(icon("table"), "Protein expression table"),
    solidHeader = TRUE,
    status = "primary",
    width = NULL,
    uiOutput("emptyTable")
  ),
  tabBox(
    title = "",
    width = NULL,
    tabPanel(
      title = tagList(icon("bar-chart"), "Intensity Distribution"),
      uiOutput("sampleDistributionBoxPanel")
    ),
    tabPanel(
      title = tagList(icon("area-chart"), "Density Plot"),
      uiOutput("sampleDistributionDensityPanel")
    ),
    tabPanel(title = tagList(icon("object-group"), "PCA"),
             uiOutput("pcaUI")),
    tabPanel(title = tagList(icon("sitemap"), "Sample Correlation"),
             uiOutput("dendUI"))
  )
)))
