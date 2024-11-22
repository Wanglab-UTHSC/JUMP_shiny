# ui-data-import.R

fluidPage(
  includeCSS("www/style/theme.css"),
  fluidRow(column(
  3,
  tabBox(
    title = "",
    id = "datasource",
    width = NULL,
    tabPanel(
      tagList(icon("cloud-upload"), "Upload"),
      radioButtons(
        "dataType",
        "Choose your file type",
        choices = c("jumpshiny",
                    "jumpq",
                    "jump_batch")
      ),
      fileInput(
        "uploadExpressionData",
        "Upload Raw Data",
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
    tabPanel(
      tagList(icon("folder-open"), "Example"),
      uiOutput("dataSourceSelect"),
      tags$p("Example JUMP raw data for illustration, containing 11 signals."),
      
      do.call(actionBttn, c(
        list(
          inputId = "ExpressionDataSample",
          label = "Import data",
          icon = icon("play")
        ),
        actionBttnParams
      ))
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
    a(href = "data/meta_DIA.zip",
      "Download example sample information table",
      download = NA,
      target = "_blank"
    ),
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
    title = tagList(icon("info-circle"), "Summary"),
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
      title = tagList(icon("bar-chart"), "Intensity distribution"),
      uiOutput("sampleDistributionBoxPanel")
    ),
    # tabPanel(
    #   title = tagList(icon("filter"), "Filtering Threshold"),
    #   uiOutput("lowCountFilterByCutoffUI")
    # ),
    tabPanel(
      title = tagList(icon("area-chart"), "Density plot"),
      uiOutput("sampleDistributionDensityPanel")
    ),
    tabPanel(title = tagList(icon("object-group"), "QQ plot"),
             uiOutput("QQplot")),
    tabPanel(title = tagList(icon("object-group"), "PCA"),
             uiOutput("pcaUI")),
    tabPanel(title = tagList(icon("sitemap"), "Sample-sample distance"),
             uiOutput("dendUI"))
  )
)))
