#server-enrichment.R

library(clusterProfiler)
library(R.utils)
library(enrichplot)
library(ggnewscale)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
library(plotly)
library(data.table)
library(DT)


options(clusterProfiler.download.method = "wget")

FilterRun <- reactiveValues(FilterRunValue = FALSE)
enRun <- reactiveValues(enRunValue = FALSE)
customPW <- reactiveValues(term2gene = NULL, term2name = NULL)


.clean_genes <- function(x) {
  x <- trimws(x); x <- x[nzchar(x)]
  toupper(x)
}
.try_read_maybe_tab_or_csv <- function(p) {
  dt <- tryCatch(fread(p, sep="\t", header=TRUE, data.table=FALSE), error=function(e) NULL)
  if (is.null(dt) || ncol(dt) == 1) {
    dt <- fread(p, sep=",", header=TRUE, data.table=FALSE)
  }
  dt
}

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
          "P-val" = "en_top"
        ),
        justified = TRUE,
        status = "primary"
      ),
      uiOutput("enrich_sigLevel"),
      numericInput("en_log2fcUp", 
                   label = "Upregulated Log2-fold cutoff", 
                   min = 0, 
                   max = 10, 
                   step = 0.1,
                   value = 0.3),
      numericInput("en_log2fcDown", 
                   label = "Downregulated Log2-fold cutoff", 
                   min = -10, 
                   max = 0, 
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
    "en_FDR" = tagList(numericInput(
      "enrichmentFDR",
      "FDR Cut-off",
      min = 0.01,
      max = 1,
      value = 0.01,
      step = 0.01
    ),
    textOutput("enrichmentProteinPreview")),
    "en_top" = tagList(numericInput(
      inputId = "en_topProt",
      label = "P-value cutoff",
      value = 0.05,
      min = 0,
      max = 1,
      step = 0.01
    ),
    textOutput("enrichmentProteinPreview")
  ))
})


# Preview proteins count -----
observeEvent(input$enrichmentFDR, {
  res <- variables$result
  prot_count <-
    nrow(res[res$FDR <= input$enrichmentFDR, ])
  output$enrichmentProteinPreview <- renderText({
    paste0(
      "Protein number: ",
      prot_count
    )
  })
})
observeEvent(input$en_topProt, {
  res <- variables$result
  prot_count <-
    nrow(res[res$"p-value" <= input$en_topProt, ])
  output$enrichmentProteinPreview <- renderText({
    paste0(
      "Protein number: ",
      prot_count
    )
  })
})


observeEvent(input$enrich_sigFilter,{
  
  # progressSweetAlert(
  #   session = session,
  #   id = "filterProgress",
  #   title = "Load in parameter",
  #   display_pct = TRUE,
  #   value = 0
  # )
  
  # updateProgressBar(
  #   session = session,
  #   id = "filterProgress",
  #   title = "Processing data",
  #   value = 30
  # )
  #Obtain filter parameters
  data <- variables$result
  res <- variables$res
  
  if(input$filterEN == "No"){
    sampleData = data
  }else{
    if(input$enrich_filterSelectType == "en_list"){
      selectList <- row.names(data) %in% unlist(strsplit(x = input$enrichmentTextList, split = "[\r\n]"))
      sampleData <- data[selectList,]
      
      #updateTextAreaInput(session = session, inputId = "enrichmentTextList",value = input$heatmapTextList)
    }
    else if(input$enrich_filterSelectType == "en_FDR"){
      selectList <-
        row.names(data) %in% row.names(res[res$FDR <= input$enrichmentFDR, ])
      sampleData <- data[selectList,]
      
      #updateSliderInput(session = session,inputId = "enrichmentFDR",value = input$enrichmentFDR)
    }
    else if(input$enrich_filterSelectType == "en_top"){
      pvalCut <- input$en_topProt
      # selectList <- 
      #   row.names(data) %in% row.names(res[res$`p-value`<=pvalCut,])
      sampleData <- data[which(data$'p-value' <= pvalCut),]
      
      #updateNumericInput(session = session,inputId = "en_topProt",value = input$en_topProt)
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
    
    # updateProgressBar(
    #   session = session,
    #   id = "filterProgress",
    #   title = "Data Selected",
    #   value = 50
    # )
    
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
  
  
  # 
  # updateProgressBar(
  #   session = session,
  #   id = "filterProgress",
  #   title = "Log2 fold change applied",
  #   value = 50
  # )
  
  
  variables$enrichedDataList <- sampleData
  #write.csv(sampleData,"test/filteredResult.csv")
  
  # updateProgressBar(
  #   session = session,
  #   id = "filterProgress",
  #   title = "Output data table",
  #   value = 80
  # )
  
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
        selection = 'single',
        extensions = c("Scroller", "Buttons"),
        option = list(
          dom = 'lfrtip',
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
  
  
  # updateProgressBar(
  #   session = session,
  #   id = "filterProgress",
  #   title = "Data Saved",
  #   value = 100
  # )
  
  FilterRun$FilterRunValue <- input$enrich_sigFilter
  
  # closeSweetAlert(session = session)
  # sendSweetAlert(
  #   session = session,
  #   title = "DONE",
  #   text = "Data filtering done.",
  #   type = "success"
  # )
  
})



output$preEnrichResultTable <- renderUI({
  if(FilterRun$FilterRunValue){
    tagList(fluidRow(column(
      12,
      downloadButton("downLoadFilTable", "Download Filtered Table (CSV)"),
      tags$br(),
      DT::dataTableOutput('preEnrichmentTable') %>% withSpinner()
    ))
    )} else {
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
observeEvent(input$useCustomGenes,{
  data <- variables$result
  
  customGenes <- unlist(strsplit(x = input$customGenes, split = "[\r\n]"))
  
  variables$customGenes <- customGenes
  matched <- intersect(customGenes,data$GN)
  genes_notfound <- setdiff(customGenes,data$GN)
  sendSweetAlert(
    session = session,
    title = "Success!",
    text = "Custom gene list uploaded.",
    type = "success"
  )
  showNotification("Custom gene list uploaded", type = "message", duration = 10)
  
})


#Helper funciton to identify the symbols user entered
# Normalize common ID quirks (version suffixes, isoforms)
normalize_ids <- function(ids) {
  ids <- trimws(ids)
  ids <- gsub("\\.\\d+$", "", ids)   # drop version: ENSG..., NM_... .1
  ids <- toupper(ids)               # UniProt often uppercase
  unique(ids[nzchar(ids)])
}
key_type <- function(customList) {
  ids <- customList
  
  # Patterns:
  # Ensembl genes: human(most): ENSG..., mouse: ENSMUSG..., general: ENS[A-Z]{0,3}G\d+
  is_ensembl <- grepl("^ENS[A-Z]{0,3}G\\d+$", ids)
  
  # RefSeq: NM_, NR_, NP_ (curated), XM_, XR_, XP_ (predicted), NC_ (genomic)
  is_refseq  <- grepl("^(NM|NR|NP|XM|XR|XP|NC)_[0-9]+(\\.[0-9]+)?$", ids)
  
  # UniProt accessions:
  #  Classic 6-char: [OPQ][0-9][A-Z0-9]{3}[0-9]  or  [A-NR-Z][0-9][A-Z0-9]{3}[0-9]
  #  New style: A0Axxxxx (A0A + 7 alnum), allow -isoform suffixes
  is_uniprot <- grepl("^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|A0A[A-Z0-9]{7})(-[0-9]+)?$", ids)
  
  # Entrez ID
  is_entrez  <- grepl("^[0-9]+$", ids)
  
  # Compute match rates
  n <- length(ids)
  rate <- c(
    ENSEMBL = sum(is_ensembl)/n,
    REFSEQ  = sum(is_refseq)/n,
    UNIPROT = sum(is_uniprot)/n,
    ENTREZID= sum(is_entrez)/n
  )

  
  best <- names(rate)[which.max(rate)]
  confidence <- max(rate)
  if (confidence == 0) {
    best <- "SYMBOL"
    return(best)
  } else {
    return(best = best)
  }
}

observeEvent(input$customPathwayFile, {
  #req(input$customPathwayFile)
  path <- input$customPathwayFile$datapath
  fmt  <- input$pathwayFileFormat %||% "std"
  
  term2gene <- NULL; term2name <- NULL
  
  if (fmt == "std") {
    dt <- .try_read_maybe_tab_or_csv(path)
    shiny::validate(need(!is.null(dt), "Could not read file."))
    shiny::validate(need(all(c("Pathway","Description","Gene") %in% colnames(dt)),
                  "Standard format needs columns: Pathway, Description, Gene."))
    dt$Gene <- .clean_genes(dt$Gene)
    term2gene <- unique(dt[, c("Pathway","Gene")]); colnames(term2gene) <- c("TERM","GENE")
    term2name <- unique(dt[, c("Pathway","Description")]); colnames(term2name) <- c("TERM","NAME")
    
  } else if (fmt == "two_column") {
    dt <- .try_read_maybe_tab_or_csv(path)
    shiny::validate(need(!is.null(dt), "Could not read file."))
    if (!all(c("Pathway","Gene") %in% colnames(dt)) && ncol(dt) >= 2) {
      colnames(dt)[1:2] <- c("Pathway","Gene")
    }
    shiny::validate(need(all(c("Pathway","Gene") %in% colnames(dt)),
                  "Two-column format needs columns: Pathway, Gene."))
    dt$Gene <- .clean_genes(dt$Gene)
    term2gene <- unique(dt[, c("Pathway","Gene")]); colnames(term2gene) <- c("TERM","GENE")
    
  } else if (fmt == "multi_column") {
    dt <- .try_read_maybe_tab_or_csv(path)
    colnames(dt)[1] <- "Pathway"
    long <- data.frame(Pathway = rep(dt$Pathway, ncol(dt)-1),
                       Gene = as.vector(as.matrix(dt[, -1, drop=FALSE])),
                       stringsAsFactors = FALSE)
    long$Gene <- .clean_genes(long$Gene)
    term2gene <- unique(subset(long, nzchar(Gene))[, c("Pathway","Gene")])
    colnames(term2gene) <- c("TERM","GENE")
  }
  
  #shiny::validate(need(!is.null(term2gene) && nrow(term2gene) > 0, "No valid TERM–GENE pairs found."))
  
  term2gene <- subset(unique(term2gene), nzchar(TERM) & nzchar(GENE))
  customPW$term2gene <- term2gene
  customPW$term2name <- term2name
  
  # Quiet success ping (no details leaked)
  updateSelectInput(session, "method", selected = "custom_enrich")
  showNotification("Switched enrcihment method to 'Custom Annotation(uploaded)'.",
                   type = "message", duration = 4)
  showNotification("Custom pathways loaded.", type = "message", duration = 4)
})
observeEvent(input$runEnrichmentAnalysis,{
  
  #Load Data
  data <- variables$enrichedDataList
  wholeData <- variables$CountData
  geneLst <- wholeData$GN
  if(length(variables$customGenes)!=0){
    gene <- variables$customGenes
    key = key_type(gene)
    
  }else{
    gene <- data$GN
    key = key_type(gene)
  }
  #method selection
  method <- input$method
  
  #ontology
  ont <- input$ont
  
  #obtain other parameters
  org = input$orgSelect
  
  
  #Obtain whether key is ENSEMBL, GN, or PID
  
  pval = input$enrich_pvalCut
  qval = input$enrich_qvalCut
  minGSSize <- input$minGSSize
  maxGSSize <- input$maxGSSize
  
  # updateProgressBar(
  #   session = session,
  #   id = "enrichmentProgress",
  #   title = "Convert keytypes",
  #   value = 50
  # )
  
  #translate to symbols
  

  
  
  #translate kegg code
  kegg_org = ""
  if(org == "org.Mm.eg.db"){
    kegg_org = "mmu"
  }else if(org == "org.Hs.eg.db"){
    kegg_org = "hsa"
  }else{
    kegg_org = "rno"
  }
  
  gene_go <- gene[!duplicated(gene)]
  geneLst_go <- geneLst[!duplicated(geneLst)]
  
  showNotification("Enrichment analysis starting...", type = "message",duration = 50)
  if(method == "go_enrich"){
    ego <- try(
      enrichGO(gene          = gene_go,
               universe      = geneLst_go,
               OrgDb         = org,
               keyType       = key,
               ont           = ont,
               pvalueCutoff  = pval,
               qvalueCutoff  = qval,
               minGSSize     = minGSSize,
               maxGSSize     = maxGSSize,
               readable      = TRUE),
      silent = TRUE
    )

    if (inherits(ego, "try-error")) {
      sendSweetAlert(session, "Enrichment failed", "GO enrichment error. Please check organism/parameters.", type = "error")
      return()
    }
    variables$en_result <- ego
  }
  else if(method == "gsea_enrich"){
    geneList <- data$Log2Fold
    name <- row.names(data)
    geneList2 <- as.data.frame(cbind(name,geneList))
    if(key != "ENTREZID") name <- bitr(name, fromType=key, toType="ENTREZID", OrgDb=org)
    
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
    if(key != "ENTREZID"){
      gene_kegg <- bitr(gene, fromType=key, toType="ENTREZID", OrgDb=org)
    }
    
    
    ego <- enrichKEGG(gene  = gene_kegg$ENTREZID,
                      organism = kegg_org,
                      minGSSize = minGSSize,
                      maxGSSize = maxGSSize,
                      pvalueCutoff  = pval,
                      qvalueCutoff  = qval)
    ego_readable <- setReadable(ego,OrgDb = org,keyType = "ENTREZID")
    
    variables$en_result <- ego_readable
  }
  else if (method == "custom_enrich"){
    term2gene <- customPW$term2gene
    shiny::validate(need(!is.null(term2gene) && nrow(term2gene) > 0,
                  "Please upload a custom pathway database first."))
    
    # Intersect user genes to uploaded gene universe (IDs must match as text)
    gene_use <- intersect(.clean_genes(gene), unique(term2gene$GENE))
    shiny::validate(need(length(gene_use) > 0,
                  "None of your genes matched the uploaded pathways. Check ID types/case."))
    
    # Optional background universe from your dataset (keeps p-values conservative)
    bg <- intersect(.clean_genes(variables$CountData$GN), unique(term2gene$GENE))
    if (length(bg) < length(gene_use)) bg <- unique(term2gene$GENE)  # fallback to all uploaded genes
    
    ego <- enricher(
      gene         = gene_go,
      TERM2GENE    = term2gene,
      TERM2NAME    = customPW$term2name,   # can be NULL
      pvalueCutoff = pval,
      qvalueCutoff = qval,
      minGSSize    = minGSSize,
      maxGSSize    = maxGSSize,
      universe     = bg
    )
    variables$en_result <- ego
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
  
  
  # updateProgressBar(
  #   session = session,
  #   id = "enrichmentProgress",
  #   title = "Output Result",
  #   value = 90
  # )
  
  shinyCatch(
    variables$en_table <- variables$en_result@result,
    blocking_level = "error"
    )
  
  output$downloadEnResultTable <- downloadHandler(
    filename = "enrichment_result.csv",
    content = function(file) {
      write.csv(variables$en_table,file)
    }
  )
  
  table <- variables$en_table
  #table <- fortify(ego,showCategory = length(ego@result$ID))
  output$afterEnrichResultTable <- DT::renderDataTable({ 
    data <- table
    data$`pvalue` = formatC(data$`pvalue`, format = "e", digits = 3)
    data$`qvalue` = formatC(data$`qvalue`, format = "e", digits = 3)
    data$`p.adjust` = formatC(data$`p.adjust`, format = "e", digits = 3)
    DT::datatable(data, 
                  options = list( 
                    dom = 'lfrtip', 
                    scrollX = TRUE,
                    scrollY = 400,
                    scroller = TRUE,
                    pageLength = 5, 
                    searchHighlight = TRUE, 
                    autoWidth = TRUE 
                  ), 
                  rownames = FALSE 
    ) 
  },server = TRUE)
  
  enRun$enRunValue <- TRUE
  
  # updateProgressBar(
  #   session = session,
  #   id = "enrichmentProgress",
  #   title = "Output Result",
  #   value = 100
  # )
  # 
  # closeSweetAlert(session = session)
  # sendSweetAlert(
  #   session = session,
  #   title = "DONE",
  #   text = "Enrichment analysis done.",
  #   type = "success"
  # )
  
})



output$afterEnrichResultUI <- renderUI({
  if(enRun$enRunValue){
    tagList(fluidRow(column(
      12,
      downloadButton("downloadEnResultTable", "Download Enrichment Result Table (CSV)"),
      tags$br(),
      DT::dataTableOutput('afterEnrichResultTable') %>% withSpinner()
    )))
  }else{
    helpText("Click [Run Enrichment Analysis] to obtain enrichment results.")
  }
})


#-------------Plot Data---------------
#Bubble Plot

output$bubbleplotObject <- renderPlot({
  shiny::validate(
    shiny::need(nrow(variables$en_result)!=0, message = "No enriched term found")
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
      plotOutput("bubbleplotObject",height = 800)%>% withSpinner()
    ))
  }else if(enRun$enRunValue && nrow(variables$en_result) == 0){
    helpText("No enriched term found.")
  }else{
    helpText("No data to plot. Run enrichment analysis first.")
  }
})



output$enrichmentMap <- renderPlot({
  shiny::validate(
    shiny::need(nrow(variables$en_result)!=0, message = "No enriched term found")
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
      plotOutput("enrichmentMap",height = 800)%>% withSpinner()
    ))
  }else if(enRun$enRunValue && nrow(variables$en_result) == 0){
    helpText("No enriched term found.")
  }
  else{
    helpText("No data to plot. Run enrichment analysis first.")
  }
})
