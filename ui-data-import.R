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
        "uploadCountData",
        "Upload Raw Data",
        accept = c("text/csv",
                   "text/comma-separated-values,text/plain",
                   ".csv"),
        buttonLabel = "Upload...",
        placeholder = "No file has been uploaded."
      ),
      helpText("Text file in .tsv/.csv format, and the first column should be Accession Number, second column should be Gene Names, Third should be Description.")
    ),
    tabPanel(
      tagList(icon("folder-open"), "Example"),
      uiOutput("dataSourceSelect"),
      tags$p("Example JUMP raw data for illustration, containing 11 signals."),
      
      do.call(actionBttn, c(
        list(
          inputId = "CountDataSample",
          label = "1. Import Data",
          icon = icon("play")
        ),
        actionBttnParams
      ))
    )
  ),
  box(
    title = tagList(icon("tags"), "Group Assignment"),
    solidHeader = TRUE,
    status = "primary",
    width = NULL,
    fileInput(
      "uploadGroup",
      "Upload Group Information",
      accept = c("text/csv",
                 "text/comma-separated-values,text/plain",
                 ".csv"),
      buttonLabel = "Upload...",
      placeholder = "No file has been uploaded."
    ),
    helpText("Text file in .tsv/.csv format, and the first column should be named as Sample."),
    do.call(actionBttn, c(
      list(
        inputId = "confirmedGroupList",
        label = "2. Assign Group Label",
        icon = icon("play")),
      actionBttnParams
      )
    ),
    footer = helpText("JUMP Suite expect first label should be Group1 (G1) and the next be Group2 (G2), and so on.")
  ),
  box(
    title = tagList(icon("layer-group"), "Group Selection"),
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
    title = tagList(icon("table"), "Protein Expression Table"),
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
    # tabPanel(
    #   title = tagList(icon("filter"), "Filtering Threshold"),
    #   uiOutput("lowCountFilterByCutoffUI")
    # ),
    tabPanel(
      title = tagList(icon("area-chart"), "Density Plot"),
      uiOutput("sampleDistributionDensityPanel")
    ),
    tabPanel(title = tagList(icon("object-group"), "MDS Plot"),
             uiOutput("mdsUI")),
    tabPanel(title = tagList(icon("object-group"), "PCA"),
             uiOutput("pcaUI")),
    tabPanel(title = tagList(icon("sitemap"), "Hierarchical Clustering"),
             uiOutput("dendUI"))
  )
)))
