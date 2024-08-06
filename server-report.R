# server-report.R

library(officer)

runReport <- reactiveValues(runReportValue = FALSE)


# Click the generate report button, render report ----
observeEvent(input$generateReport, {
  
  # src <- normalizePath('Plot_Report.Rmd')
  # 
  # 
  # owd <- setwd(tempdir())
  # on.exit(setwd(owd))
  # file.copy(src, 'Plot_Report.Rmd', overwrite = TRUE)
  
  library(rmarkdown)

  
    #Summary of raw data
  rawdata <- variables$CountData
  nProteins <- nrow(rawdata)
  projN <- input$projName
  speciesN <- input$speciesN
  tissueType <- input$tissueType
  
  
  reportParameter <- list(
    nProteins = nProteins,
    
    #original data
    rawdata = variables$CountData,
    group = variables$group,
    groupList = variables$groupList,
    groupListConvert = variables$groupListConvert,
    
    batchList = variables$batchGroupList,
    batchListConvert = variables$batchGroupListConvert,
    
    #result after differential expression
    result = variables$result,
    # res = variables$res,
    
    #result of normalization
    norData = variables$normedData,
    
    #result of enrichment
    en_result = variables$en_result,
    
    #data import section
    sampleDistributionBar = variables$sampleDistributionBar,
    sampleDistributionDensity = variables$sampleDistributionDensity,
    original_heatmap = variables$original_heatmap,
    pcaParameter = variables$pcaParameter,
    screePlot = variables$screePlot,
    pca3d = variables$pca3d,
    pca2d = variables$pca2d,
    summaryPCA = variables$summaryPCA,
    
    
    #normlization result
    norSampleDistributionBar = variables$normed_sampleDistributionBar,
    norSampleDistributionDensity = variables$normed_sampleDistributionDensity,
    norData = variables$normedData,
    normed_heatmap = variables$normed_heatmap,
    normedScreePlot = variables$normedScreePlot,
    normedPCA3d = variables$normed_pca3d,
    normedPCA2d = variables$normed_pca2d,
    normed_pcaParameter = variables$normed_pcaParameter,
    
    VolcanoPlotObject = variables$volcano,
    DE_heatmapObject = variables$de_heatmap,
    
    dotplot = variables$dotplot,
    enrichmentMap = variables$enrichMap,
    
    #session info
    projName = projN,
    speciesName = speciesN,
    tissueT = tissueType,
    
    nSamples = variables$groupText,
    
    rendered_by_shiny = TRUE,
    
    #covariates correction
    aftercov = variables$corrected_data
    
    
    # mdsPlot = NA,
    # mdsPlotplot = NA
    
  )
  
  variables$params <- reportParameter


  #variables$reportFile <- out
  runReport$runReportValue <- input$generateReport

})

output$renderReportButton <- renderUI({
  if (runReport$runReportValue) {
    tagList(tags$br(), downloadButton('downloadReport',label = "Download Report"))
  } else {
    helpText("Click [Generate Report] to retrieve parameters.")
  }
})

output$downloadReport <- downloadHandler(
  filename = function() {
    paste('JUMP_Suite_Report', sep = '.', switch(
      input$format,
      HTML = 'html',
      docx = "docx"
    ))
  },
  
  content = function(file) {
  
    withProgress(message = 'Rendering report...',{
      tempReport = file.path(tempdir(), "Report.Rmd")
      file.copy("report/Report.Rmd", tempReport, overwrite = TRUE)
      render(tempReport,output_file = file, params = variables$params, switch(
        input$format,
        HTML = html_document(),
        docx = word_document()
      ),envir = new.env(parent = globalenv()))
    })

    
  }
)

output$renderDownloadButton <- renderUI({
  if (runReport$runReportValue) {
    tagList(tags$br(), downloadButton('downloadTBReport',label = "Download All Tables"))
  } else {
    helpText("Click [Generate Report] for generation.")
  }
})

output$downloadTBReport <- downloadHandler(
  # filename = function() {
  #   paste('JUMP_Suite_Report', sep = '.', switch(
  #     input$format,
  #     HTML = 'html'
  #   ))
  # },
  filename = function(){
    paste("JUMP_Suite_", Sys.Date(), ".zip", sep = "")
  },
  content = function(file){
    reportVar <- reactiveValues(
      rawdata = variables$CountData,
      diffExpr_result = variables$result,
      normalizedData = variables$normedData,
      enrichment_result = variables$en_table,
      covariance_correction = variables$corrected_data
    )
    
    temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
    dir.create(temp_directory)
    
    reactiveValuesToList(reportVar) %>%
      imap(function(x,y){
        if(!is.null(x)){
          file_name <- glue("{y}.csv")
          readr::write_csv(x, file.path(temp_directory, file_name))
        }
      })
    
    
    zip::zip(
      zipfile = file,
      files = dir(temp_directory),
      root = temp_directory
    )
    
    
    
  },
  contentType = "application/zip"
)



# # Tab click logs -----
# observeEvent(input$sider, {
#   clickTab <- switch(
#     input$sider,
#     "welcome" = "Guidance",
#     "dateImport" = "Data Import",
#     "calculationTab" = "Calculation",
#     "maplotTab" = "MA Plot",
#     "volcanoplotTab" = "Volcano Plot",
#     "pcaTab" = "PCA Analysis",
#     "heatmapTab" = "Heatmap",
#     "expressionTab" = "Expression",
#     "reportTab" = "Report"
#   )
#   variables$logList <-
#     rbind(variables$logList,
#           list(
#             "Time" = as.character(Sys.time()),
#             "Type" = "Tab",
#             "Action" = paste0(clickTab, collapse = "-"),
#             "Parameters" = ""
#           ),
#           stringsAsFactors = FALSE)
# })
# 
# 
# # Load Sample Data botton log -----
# observeEvent(input$CountDataSample, {
#   variables$logList <- rbind(
#     variables$logList,
#     list(
#       "Time" = as.character(Sys.time()),
#       "Type" = "Button",
#       "Action" = "Load sample data",
#       "Parameters" = input$SampleDatabase
#     ),
#     stringsAsFactors = FALSE
#   )
# })
# 
# # Click TCC botton log ----
# observeEvent(input$TCC, {
#   TCCParaLog <- paste(
#     "Filter low count genes threshold:",
#     input$filterLowCount,
#     "Normalization method:",
#     input$normMethod,
#     "DEGs identify method:",
#     input$testMethod,
#     "Interation:",
#     input$iteration,
#     "FDR:",
#     input$fdr,
#     "Elimination of Potential DEGs:",
#     input$floorpdeg,
#     sep = " "
#   )
#   variables$logList <- rbind(
#     variables$logList,
#     list(
#       "Time" = as.character(Sys.time()),
#       "Type" = "Button",
#       "Action" = "Run TCC",
#       "Parameters" = TCCParaLog
#     ),
#     stringsAsFactors = FALSE
#   )
# })
# 
# 
# # Click MA botton log ----
# observeEvent(input$makeMAPlot , {
#   MAParaLog <- paste(
#     "Point Size:",
#     input$pointSize,
#     "FDR:",
#     input$maFDR,
#     "DEGs color:",
#     input$fdrColor,
#     sep = " "
#   )
#   variables$logList <- rbind(
#     variables$logList,
#     list(
#       "Time" = as.character(Sys.time()),
#       "Type" = "Button",
#       "Action" = "Generate MA-Plot",
#       "Parameters" = MAParaLog
#     ),
#     stringsAsFactors = FALSE
#   )
# })
# 
# # Click Volcano botton log ----
# 
# observeEvent(input$makeVolcanoPlot , {
#   VolcanoParaLog <- paste(
#     "Fold Change cut-off:",
#     paste(input$CutFC, collapse = "~"),
#     "p-value cut-off:",
#     input$Cutpvalue,
#     "Point Size:",
#     input$pointSize,
#     "Down-regulate:",
#     input$downColor,
#     "Up-regulate:",
#     input$upColor,
#     sep = " "
#   )
#   variables$logList <- rbind(
#     variables$logList,
#     list(
#       "Time" = as.character(Sys.time()),
#       "Type" = "Button",
#       "Action" = "Generate Volcano Plot",
#       "Parameters" = VolcanoParaLog
#     ),
#     stringsAsFactors = FALSE
#   )
# })
# 
# 
# # Click PCA botton log ----
# 
# observeEvent(input$pcRun, {
#   if(input$pcFDR != ""){
#     pcaFDR <- paste0("FDR:", input$pcFDR, sep = " ")
#   } else {
#     pcaFDR <- ""
#   }
#   pcaParaLog <- paste(
#     pcaFDR,
#     "Center:",
#     input$pcCenter,
#     "Scale:",
#     input$pcScale,
#     "Log transform:",
#     input$pcTransform,
#     "Source:",
#     input$pcData,
#     "Hierarchical Clustering Method:",
#     input$dendMethod,
#     sep = " "
#   )
#   variables$logList <- rbind(
#     variables$logList,
#     list(
#       "Time" = as.character(Sys.time()),
#       "Type" = "Button",
#       "Action" = "Run PCA Analysis",
#       "Parameters" = pcaParaLog
#     ),
#     stringsAsFactors = FALSE
#   )
# })

# # Input table ----
# output$inputLogTable <- DT::renderDataTable({
#   DT::datatable(variables$logList)
# })

# AllInputs <- reactive({
#   x <- reactiveValuesToList(input)
#   data.frame(
#     names = names(x),
#     values = unlist(x, use.names = FALSE)
#   )
# })
# 
# output$showInputs <- renderTable({
#   AllInputs()
# })