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
        textAreaInput(
          "batchSelectViaText",
          "Input your different batch info",
          rows = 6,
          placeholder = paste0(
            "Please input batch information at here. ",
            "--Linear Example--",
            "Sample,Batch\n",
            "sig126,Batch1",
            "sig127N,Batch1",
            "sig127C,Batch1",
            "sig128N,Batch1",
            "sig128C,Batch2",
            "sig129N,Batch2",
            "sig129C,Batch2",
            "sig130N,Batch2",
            "--Internal Example--",
            "Sample,Batch,Info",
            "sig126,Batch1,internal",
            "sig127N,Batch1",
            "sig127C,Batch1",
            "sig128N,Batch1",
            "sig128C,Batch2,internal",
            "sig129N,Batch2",
            "sig129C,Batch2",
            "sig130N,Batch2"
            # sep = "\n"
          )
        ),
        do.call(actionBttn, c(
          list(
            inputId = "confirmedBatchList",
            label = "Run Batch Normalization",
            icon = icon("play")),
          actionBttnParams
        )
        ),
        footer = helpText("JUMP expect first column to be Sample and the second one be Batch.")
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
             tabPanel(title = tagList(icon("sitemap"), "Hierarchical Clustering"),
                      uiOutput("norm_dendUI"))
           ))
  )
)
