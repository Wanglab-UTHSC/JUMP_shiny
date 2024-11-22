#ui-normalization.R
fluidPage(
  includeCSS("www/style/theme.css"),
  fluidRow(
    #create side board
    column(
      3,
      box(
        title = tagList(icon("cogs"), "Normalization Method"),
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        tagList(
          selectInput(
            inputId = "normalization",
            label = "Normalization Method",
            #select between two normalization methods
            choices = c("Internal" = "internalNM",
                        "Linear" = "linearNM")
            
          )
        )
      ),
      box(
        #User will input their batch information
        title = tagList(icon("tags"),"Batch Information"),
        tags$p("Please specify the internal information as in example in the batch info if selected."),
        solidHeader = T,
        status = "primary",
        width = NULL,
        uiOutput("batchSelect"),
        uiOutput("internalReference"),
        do.call(actionBttn, c(
          list(
            inputId = "confirmedBatchList",
            label = "Run Batch Normalization",
            icon = icon("play")),
          actionBttnParams
        )
        ),
        footer = helpText("Internal reference column should specify reference sample as 'internal'.")
      )
    ),
    column(9,
           box(
             #show the data table after normalization
             title = tagList(icon("table"), "Data Table"),
             width = NULL,
             solidHeader = TRUE,
             status = "primary",
             uiOutput("normalizationResultTable")
           ),
           tabBox(
             title = "",
             width = NULL,
             tabPanel(
               title = tagList(icon("bar-chart"), "Intensity Distribution"),
               uiOutput("norm_DistributionBoxPanel")
             ),
             tabPanel(
               title = tagList(icon("area-chart"), "Density Plot"),
               uiOutput("norm_sampleDistributionDensityPanel")
             ),
             tabPanel(title = tagList(icon("object-group"), "PCA"),
                      uiOutput("norm_pcaUI")),
             tabPanel(title = tagList(icon("sitemap"), "Sample-sample distance"),
                      uiOutput("norm_dendUI"))
           ))
  )
)
