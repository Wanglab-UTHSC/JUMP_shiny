#server-enrichment.R

library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
library(plotly)

FilterRun <- reactiveValues(FilterRunValue = FALSE)
enRun <- reactiveValues(enRunValue = FALSE)

#UI of the first three tabs
output$enrich_filter <- renderUI({
  if(input$filterEN == "Yes"){
    tagList(
      radioGroupButtons(
        inputId = "enrich_filterSelectType",
        label = "Significance selection",
        choices = c(
          "List" = "en_list",
          # "By name",
          "FDR" = "en_FDR",
          "p-val" = "en_top"
        ),
        justified = TRUE,
        status = "primary"
      ),
      uiOutput("enrich_sigLevel"),
      numericInput("en_log2fcUp", 
                   label = "Upregulated Log2-fold cutoff", 
                   min = -2, 
                   max = 2, 
                   step = 0.1,
                   value = 0.3),
      numericInput("en_log2fcDown", 
                   label = "Downregulated Log2-fold cutoff", 
                   min = -2, 
                   max = 2, 
                   step = 0.1,
                   value = -0.3)
    )

  }
})

output$enrich_sigLevel <- renderUI({
  data <- variables$res
  nprotein <- nrow(data)
  
  switch(
    input$enrich_filterSelectType,
    "en_list" = textAreaInput(
      "enrichmentTextList",
      "Paste Protein List",
      rows = 5,
      placeholder = "Input protein's name (first column in the dataset), one protein id per line."
    ),
    "en_FDR" = tagList(
      tagList(
        sliderInput(
          "enrichmentFDR",
          "FDR Cut-off",
          min = 0.01,
          max = 1,
          value = 0.01
        ),
        textOutput("enrichmentProteinPreview")
      )
    ),
    "en_top" = numericInput(
      inputId = "en_topProt",
      label = "P-value cutoff",
      value = 0.05,
      min = 0,
      max = 1,
      step = 0.01
    )
  )
})


# Preview proteins count -----
observeEvent(input$enrichmentFDR, {
  res <- variables$result
  prot_count <-
    nrow(res[res$FDR <= input$enrichmentFDR, ])
  output$enrichmentProteinPreview <- renderText({
    paste0(
      "Protein number: ",
      prot_count,
      " | Generation time: ~",
      round(prot_count / 30, 2),
      "s"
    )
  })
})


observeEvent(input$enrich_sigFilter,{
  
  progressSweetAlert(
    session = session,
    id = "filterProgress",
    title = "Load in parameter",
    display_pct = TRUE,
    value = 0
  )
  
  updateProgressBar(
    session = session,
    id = "filterProgress",
    title = "Processing data",
    value = 30
  )
  #Obtain filter parameters
  data <- variables$result
  res <- variables$res
  
  if(input$filterEN == "No"){
    sampleData = data
  }else{
    if(input$enrich_filterSelectType == "en_list"){
      selectList <- row.names(data) %in% unlist(strsplit(x = input$enrichmentTextList, split = "[\r\n]"))
      sampleData <- data[selectList,]
      
      updateTextAreaInput(session = session, inputId = "enrichmentTextList",value = input$heatmapTextList)
    }
    else if(input$enrich_filterSelectType == "en_FDR"){
      selectList <-
        row.names(data) %in% row.names(res[res$FDR <= input$enrichmentFDR, ])
      sampleData <- data[selectList,]
      
      updateSliderInput(session = session,inputId = "enrichmentFDR",value = input$enrichmentFDR)
    }
    else if(input$enrich_filterSelectType == "en_top"){
      pvalCut <- input$en_topProt
      # selectList <- 
      #   row.names(data) %in% row.names(res[res$`p-value`<=pvalCut,])
      sampleData <- data[which(data$'p-value' <= pvalCut),]
      
      updateNumericInput(session = session,inputId = "en_topProt",value = input$en_topProt)
    }
    else{
      sendSweetAlert(
        session = session,
        title = "ERROR",
        text = "Select a method to filter data first!",
        type = "info"
      )
      return()
    }
    
    updateProgressBar(
      session = session,
      id = "filterProgress",
      title = "Data Selected",
      value = 50
    )
    
    upCut <- input$en_log2fcUp
    downCut <- input$en_log2fcDown
    
    tryCatch(
      {
        DownSample <- sampleData[which(sampleData$Log2Fold <= downCut),]
        UpSample <- sampleData[which(sampleData$Log2Fold >= upCut),]
        showNotification("Log2-Fold Change cut applied", type = "message")
      },
      error = function(e) {
        sendSweetAlert(
          session = session,
          title = "Input Log2FC cut out of boundary",
          text = as.character(message(e)),
          type = "error"
        )
        return()
      },
      warning = function(w) {
        sendSweetAlert(
          session = session,
          title = "Input Log2FC cut boundary warning!",
          text = "The input cutoff might be out of boundary. Please check your cutoff.",
          type = "warning"
        )
        return()
      }
    )
    
    sampleData <- rbind(DownSample,UpSample)
  }
  
  
  
  updateProgressBar(
    session = session,
    id = "filterProgress",
    title = "Log2 fold change applied",
    value = 50
  )
  
  
  variables$enrichedDataList <- sampleData
  #write.csv(sampleData,"test/filteredResult.csv")
  
  updateProgressBar(
    session = session,
    id = "filterProgress",
    title = "Output data table",
    value = 80
  )
  
  output$downLoadFilTable <- downloadHandler(
    filename = "pre_enrichment_filter.csv",
    content = function(file) {
      write.csv(sampleData, file)
    }
  )
  output$preEnrichmentTable <- DT::renderDataTable({
    if (nrow(variables$enrichedDataList) == 0) {
      DT::datatable(variables$enrichedDataList)
    } else{
      data <- variables$enrichedDataList
      colInd = grep('^sig[0-9]{3}', colnames(data))
      #data[, colInd] = round(data[, colInd], digits = 2)
      data$`p-value` = formatC(data$`p-value`, format = "e", digits = 3)
      data$FDR = formatC(data$FDR, format = "e", digits = 3)
      log2Ind = grep("Log2Fold", colnames(data))
      data[, log2Ind] = round(data[, log2Ind], digits = 3)
      
      DT::datatable(
        data,
        filter = "bottom",
        colnames = c("Uniprot ID" = 1),
        # caption = tags$caption(
        #   tags$li("Filter proteins by typing condictions (such as 2...5) in the filter boxes to filter numeric columns. ",
        #           tags$b("Copy"),
        #           ", ",
        #           tags$b("Print"),
        #           " and ",
        #           tags$b("Download"),
        #           " the filtered result for further analysis."
        #   ),
        #   tags$li(
        #     HTML("<font color=\"#B22222\"><b>Protein ID</b></font> is colored according to FDR cut-off.")
        #   )
        # ),
        selection = 'single',
        extensions = c("Scroller", "Buttons"),
        option = list(
          dom = 'lfrtip',
          # buttons =
          #   list(
          #     'copy',
          #     'print',
          #     list(
          #       extend = 'collection',
          #       buttons = c('csv', 'excel', 'pdf'),
          #       text = 'Download'
          #     )
          #   ),
          deferRender = TRUE,
          scrollY = 400,
          scrollX = TRUE,
          scroller = TRUE,
          
          pageLength = 5,
          searchHighlight = TRUE,
          orderClasses = TRUE,
          columnDefs = list(
            list(visible = TRUE, targets = -1)
          )
        )
      )  
    }
  },server = TRUE)
  
  
  updateProgressBar(
    session = session,
    id = "filterProgress",
    title = "Data Saved",
    value = 100
  )
  
  FilterRun$FilterRunValue <- input$enrich_sigFilter
  
  closeSweetAlert(session = session)
  sendSweetAlert(
    session = session,
    title = "DONE",
    text = "Data filtering done.",
    type = "success"
  )
  
})



output$preEnrichResultTable <- renderUI({
  if(FilterRun$FilterRunValue){
    tagList(fluidRow(column(
      12,
      downloadButton("downLoadFilTable", "Download Filtered Table (CSV)")
    )),
    tags$br(),
    fluidRow(column(
      12, 
      DT::dataTableOutput('preEnrichmentTable') %>% withSpinner()
    )))} else {
      helpText("Click [Run Significance Filter] to obtain Filtered Table.")
    }
})


#Parameter selection if enrich GO or GSE is selected
output$go_enrich_select <- renderUI({
  method <- input$method
  if(method == "go_enrich"||method =="gsea_enrich"){
    selectInput("ont", "Ontology",
                choices = c(
                  "ALL" = "ALL",
                  "CC" = "CC",
                  "BP" = "BP",
                  "MF" = "MF"
                ), 
                selected = "ALL")
  }
  
})



#------------Run Enrichment Test---------

observeEvent(input$runEnrichmentAnalysis,{
  
  progressSweetAlert(
    session = session,
    id = "enrichmentProgress",
    title = "Loading parameter",
    display_pct = TRUE,
    value = 0
  )
  
  updateProgressBar(
    session = session,
    id = "enrichmentProgress",
    title = "Loading data",
    value = 30
  )
  
  #Load Data
  data <- variables$enrichedDataList
  wholeData <- variables$result
  gene <- row.names(data)
  geneLst <- row.names(wholeData)
  
  #method selection
  method <- input$method
  
  #ontology
  ont <- input$ont
  
  #obtain other parameters
  org = input$orgSelect
  
  
  #Obtain whether key is ENSEMBL, GN, or PID
  key = "UNIPROT"
  pval = input$enrich_pvalCut
  qval = input$enrich_qvalCut
  minGSSize <- input$minGSSize
  maxGSSize <- input$maxGSSize
  
  updateProgressBar(
    session = session,
    id = "enrichmentProgress",
    title = "Convert keytypes",
    value = 50
  )
  
  #translate to symbols
  tryCatch(
    {
      gene_go <- bitr(gene, fromType=key, toType="SYMBOL", OrgDb=org)
      gene_go <- gene_go$SYMBOL
      gene_go <- gene_go[!duplicated(gene_go)]
      geneLst_go <- bitr(geneLst, fromType=key, toType="SYMBOL", OrgDb=org)
      geneLst_go <- geneLst_go$SYMBOL
      geneLst_go <- geneLst_go[!duplicated(geneLst_go)]
    },
    error = function(e){
      #closeSweetAlert(session = session)
      sendSweetAlert(
        session = session,
        title = "Enrichment Analysis failed",
        text = "Convert keytype failed. Please check your selected organism.",
        type = "error"
      )
      return()
    }
  )

  
  
  #translate kegg code
  kegg_org = ""
  if(org == "org.Mm.eg.db"){
    kegg_org = "mmu"
  }else if(org == "org.Hs.eg.db"){
    kegg_org = "hsa"
  }else{
    kegg_org = "rno"
  }
  
  updateProgressBar(
    session = session,
    id = "enrichmentProgress",
    title = "Perform Enrichment Analysis",
    value = 70
  )
  
  if(method == "go_enrich"){
    tryCatch({
      ego <- enrichGO(gene = gene_go,
                      universe = geneLst_go,
                      OrgDb         = org,
                      keyType       = 'SYMBOL',
                      ont           = ont,
                      pvalueCutoff  = pval,
                      qvalueCutoff  = qval,
                      readable = TRUE)
      variables$en_result <- ego
    },
      error = function(e) {
        sendSweetAlert(
          session = session,
          title = "Enrichment Analysis failed. ",
          text = "Please check your parameters.",
          type = "error"
        )
        return()
        })
  }
  else if(method == "gsea_enrich"){
    geneList <- data$Log2Fold
    name <- row.names(data)
    geneList2 <- as.data.frame(cbind(name,geneList))
    name <- bitr(name, fromType=key, toType="ENTREZID", OrgDb=org)
    geneList2 <- geneList2[geneList2$name %in% name$UNIPROT,]
    geneList <- as.numeric(geneList2$geneList)
    
    
    names(geneList) <- as.character(name$ENTREZID)
    
    geneList <- sort(geneList, decreasing = TRUE)
    
    ego <- gseGO(geneList     = geneList,
                 OrgDb        = org,
                 ont          = ont,
                 minGSSize    = minGSSize,
                 maxGSSize    = maxGSSize,
                 pvalueCutoff = pval,
                 verbose      = FALSE)
    variables$en_result <- ego
    
  }
  else if(method == "kegg_enrich"){
    gene_kegg <- bitr(gene, fromType=key, toType="ENTREZID", OrgDb=org)
    
    
    ego <- enrichKEGG(gene  = gene_kegg$ENTREZID,
                      organism = kegg_org,
                      minGSSize = minGSSize,
                      maxGSSize = maxGSSize,
                      pvalueCutoff  = pval,
                      qvalueCutoff  = qval)
    ego_readable <- setReadable(ego,OrgDb = org,keyType = "ENTREZID")
    
    variables$en_result <- ego_readable
  }
  else{
    geneList <- data$Log2Fold
    name <- row.names(data)
    geneList2 <- as.data.frame(cbind(name,geneList))
    name <- bitr(name, fromType=key, toType="ENTREZID", OrgDb=org)
    geneList2 <- geneList2[geneList2$name %in% name$UNIPROT,]
    geneList <- as.numeric(geneList2$geneList)
    
    names(geneList) <- as.character(name$ENTREZID)
    
    geneList <- sort(geneList, decreasing = TRUE)
    
    ego <- gseKEGG(geneList     = geneList,
                   organism     = kegg_org,
                   minGSSize    = minGSSize,
                   maxGSSize = maxGSSize,
                   pvalueCutoff = pval,
                   verbose      = FALSE)
    variables$en_result <- ego
    
  }
  
  
  updateProgressBar(
    session = session,
    id = "enrichmentProgress",
    title = "Output Result",
    value = 90
  )
  
  shinyCatch(
    variables$en_table <- variables$en_result@result,
    blocking_level = "error"
    )
  
  output$downloadEnResultTable <- downloadHandler(
    filename = "enrichment_result.csv",
    content = function(file) {
      write.csv(variables$en_table, file)
    }
  )
  
  output$afterEnrichResultTable <- DT::renderDataTable({
    if (nrow(variables$en_result) == 0) {
      DT::datatable(variables$en_result@result)
    }else{
      data <- variables$en_result
      data <- data@result
      data$`pvalue` = formatC(data$`pvalue`, format = "e", digits = 3)
      data$`qvalue` = formatC(data$`qvalue`, format = "e", digits = 3)
      data$`p.adjust` = formatC(data$`p.adjust`, format = "e", digits = 3)
      DT::datatable(
        data = data,
        filter = "bottom",        
        # caption = tags$caption(
        #   tags$li("Filter conditions in the bottom. ",
        #           tags$b("Copy"),
        #           ", ",
        #           tags$b("Print"),
        #           " and ",
        #           tags$b("Download"),
        #           " the enrichment result for further analysis."
        #   )
        # ),
        selection = 'single',
        extensions = c("Scroller", "Buttons"),
        option = list(
          dom = 'lfrtip',
          # buttons =
          #   list(
          #     'copy',
          #     'print',
          #     list(
          #       extend = 'collection',
          #       buttons = c('csv', 'excel', 'pdf'),
          #       text = 'Download'
          #     )
          #   ),
          deferRender = TRUE,
          scrollY = 400,
          scrollX = TRUE,
          scroller = TRUE,
          
          pageLength = 5,
          searchHighlight = TRUE,
          orderClasses = TRUE,
          columnDefs = list(
            list(visible = TRUE, targets = -1)
          )
        )
      )
    }
  },server = TRUE)
  
  enRun$enRunValue <- input$runEnrichmentAnalysis
  
  updateProgressBar(
    session = session,
    id = "enrichmentProgress",
    title = "Output Result",
    value = 100
  )
  
  closeSweetAlert(session = session)
  sendSweetAlert(
    session = session,
    title = "DONE",
    text = "Enrichment analysis done.",
    type = "success"
  )
  
})


output$afterEnrichResultUI <- renderUI({
  if(enRun$enRunValue){
    tagList(fluidRow(column(
      12,
      downloadButton("downloadEnResultTable", "Download Enrichment Result Table (CSV)")
    )),
    tags$br(),
    fluidRow(column(
      12, 
      DT::dataTableOutput('afterEnrichResultTable') %>% withSpinner()
    )))
  }else{
    helpText("Click [Run Enrichment Analysis] to obtain enrichment results.")
  }
})


#-------------Plot Data---------------
#Bubble Plot

output$bubbleplotObject <- renderPlot({
  validate(
    need(nrow(variables$en_result)!=0, message = "No enriched term found")
  )
  data <- variables$en_result
  p <- dotplot(data,showCategory = 20)+ggtitle("Dot plot")
  variables$dotplot <- p
  p
  
})

output$EN_dotPlot <- renderUI({
  if(enRun$enRunValue){
    fluidRow(column(
      12,
      plotOutput("bubbleplotObject")%>% withSpinner()
    ))
  }else if(enRun$enRunValue && nrow(variables$en_result) == 0){
    helpText("No enriched term found.")
  }else{
    helpText("No data to plot. Run enrichment analysis first.")
  }
})



output$enrichmentMap <- renderPlot({
  validate(
    need(nrow(variables$en_result)!=0, message = "No enriched term found")
  )
  data <- variables$en_result
  x <- pairwise_termsim(data)
  p <- emapplot(x)+ggtitle("Enrichment Map")
  variables$enrichMap <- p
  p
  
})

output$EN_network <- renderUI({
  if(enRun$enRunValue){
    fluidRow(column(
      12,
      plotOutput("enrichmentMap")%>% withSpinner()
    ))
  }else if(enRun$enRunValue && nrow(variables$en_result) == 0){
    helpText("No enriched term found.")
  }
  else{
    helpText("No data to plot. Run enrichment analysis first.")
  }
})
