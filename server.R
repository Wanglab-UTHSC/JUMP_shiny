library(shiny)

options(shiny.maxRequestSize = 60*1024^2)

source(file = "global.R",
       local = TRUE,
       encoding = "UTF-8")

# Define server
shinyServer(function(input, output, session) {
  
  # the modal dialog where the user can enter the query details.
  updates_modal <- modalDialog(
    title = tags$h1("240702 Updates",style="text-align:center"),
    tags$div(
      tags$h3("New Feature"),
      tags$ul(
        tags$li("Differential Expression Analysis: Now missing value imputation is available for pairwise comparison"),
        #tags$li("Heatmap under DE: Allow user to choose clustering dendrogram on rows/columns"),
        #tags$li("Covariates analysis could now work properly")
      ),
      tags$hr(),
      tags$h3("Fixes"),
      tags$ul(
        tags$li("Sample distribution boxplot filter fixed")
        #tags$li("PCA 2D plot label fixed"),
        #tags$li("Moving SD calculation approach fixed")
      ),
    ),
    easyClose = T
  )
  
  # Show the model on start up ...
  showModal(updates_modal)
  
  
  #store global variables
  variables = reactiveValues(
    dup = 0,
    simulationData = NULL,
    count.data = NULL,
    normed.count.data = NULL,
    CountData = data.frame(),
    dataToShow = data.frame(),
    groupList = NULL,
    groupListConvert = NULL,
    group_import = NULL,
    batchGroupList = NULL,
    batchGroupListConvert = NULL,
    groupText = "",
    
    #differentiation result
    result = data.frame("Results will show here." = character(0)),
    res = data.frame(),
    zeroValue = "",
    norData = "",
    runTCCCode = "",
    runMAPlot = "",
    runVolcanoPlot = "",
    runHeatmap = "",
    runPCACode = "",
    sampleDistributionBar = "",
    sampleDistributionDensity = "",
    normed_sampleDistributionDensity = "",
    normed_sampleDistributionBar = "",
    norSampleDistributionBar = "",
    norSampleDistributionDensity = "",
    MAPlotObject = "",
    VolcanoPlotObject = "",
    
    normedData = data.frame(),
    
    mdsPlot = list(),
    mdsPlotplot = NULL,
    
    data.pca = NULL,
    pcaParameter = NULL,
    normed_pcaParameter = NULL,
    screePlot = NULL,
    normedScreePlot = NULL,
    pca3d = NULL,
    pca2d = NULL,
    normed_pca3d = NULL,
    normed_pca2d = NULL,
    summaryPCA = NULL,
    
    original_heatmap = NULL,
    de_heatmap = NULL,
    
    expressionData = NULL,
    expressionLevelBar = NULL,
    expressionLevelBox = NULL,
    expressionLevelCountTable = NULL,
    expressionLevelResultTable = NULL,
    
    #KSEM
    eSEMRaw = data.frame(),
    eSEMnorm = data.frame(),
    
    #enrichment
    en_result = data.frame(),
    en_table  = data.frame(),
    
    #covariance
    corrected_data = data.frame(),
    
    logList = data.frame(
      "Time" = vector(),
      "Type" = vector(),
      "Action" = vector(),
      "Parameters" = vector()
    ),
    reportFile = NULL
  )
  source(file = "server-blockRand.R",
         local = TRUE,
         verbose = FALSE,
         encoding = "UTF-8")
  source(file = "server-data-import.R",
         local = TRUE,
         encoding = "UTF-8")
  source(file = "server-normalization.R",
         local = TRUE,
         encoding = "UTF-8")
  source(file = "server-differential-expression.R",
         local = TRUE,
         encoding = "UTF-8")
  source(file = "server-enrichment.R",
         local = TRUE,
         encoding = "UTF-8")
  source(file = "server-jumpn.R",
         local = TRUE,
         encoding = "UTF-8")
  source(file = "server-covariance.R",
         local = TRUE,
         encoding = "UTF-8")
  # source(file = "server-KSEM.R",
  #        local = TRUE,
  #        encoding = "UTF-8")
  source(file = "server-report.R",
         local = TRUE,
         encoding = "UTF-8")
  
})
