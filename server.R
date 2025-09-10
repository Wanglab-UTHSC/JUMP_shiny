#Acknowledgement: 
#Some code of this program was adapted from a program related under the MIT license:
#TCC-GUI(https://github.com/swsoyee/TCC-GUI)

library(shiny)


source(file = "global.R",
       local = TRUE,
       encoding = "UTF-8")

# Define server
shinyServer(function(input, output, session) {
  
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
    norData = "",
    runVolcanoPlot = "",
    runHeatmap = "",
    runPCACode = "",
    sampleDistributionBar = "",
    sampleDistributionDensity = "",
    normed_sampleDistributionDensity = "",
    normed_sampleDistributionBar = "",
    norSampleDistributionBar = "",
    norSampleDistributionDensity = "",
    VolcanoPlotObject = "",
    
    normedData = data.frame(),
    
    data.pca = NULL,
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

  source(file = "server-covariance.R",
         local = TRUE,
         encoding = "UTF-8")
  
})
