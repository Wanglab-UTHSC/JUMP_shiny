#server-differential-expression.R
library(scatterD3)
library(gplots)
library(ggplot2)
library(plotly)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(DT)
library(curl)
library(pheatmap)
library(gridlayout)
library(zoo)
library(shinyBS)

options(digits = 5, scipen = 10)

source(file = "statTest.R",
       local = TRUE,
       encoding = "UTF-8")
source(file = "preprocess.R",
       local = TRUE,
       encoding = "UTF-8")

#if run DE button is clicked, run the test

DERun <- reactiveValues(DERunValue = FALSE)
volrun <- reactiveValues(volrunval = FALSE)
HeatmapRun <- reactiveValues(heatmapRunValue = FALSE, height = 300)

#load group selection
metaData1 <- reactive(variables$group)
# Specification of groups of samples
output$groups2 = renderUI({
  df = metaData1()
  if (is.null(colnames(df))) {
    vars = NULL
  } else {
    vars = colnames(df)[2: ncol(df)]    
  }
  selectInput("groups2", "Grouping variable",
              choices = vars, selected = vars[1])
})

#Select pairwise/multiple compare
output$compSelect = renderUI({
  radioButtons(
    inputId = "compareSelection",
    label = "Select pairwise or multiple group comparison",
    choices = c("Pairwise" = "pair","Multiple Comparison" = "multicomp"),
    selected = "pair"
  )
})


#choose group number for comparison
output$group3 <- renderUI({
  output = tagList()
  group <- input$groups2
  df = metaData1()
  if (is.null(colnames(df)) || is.null(group)) {
    vars = NULL
  } else {
    idx = grep(group, colnames(df))
    vars = unique(df[,idx])
  }
  
  ngroup <- length(vars)
  
  if(input$compareSelection == "pair"){

      output[[1]] = selectInput(
        inputId = "group1select",
        label = "Select first group for comparison",
        choices = vars,
        selected = vars[1]
      )
      output[[2]] = selectInput(
        inputId = "group2select",
        label = "Select second group for comparison",
        choices = vars,
        selected = vars[2]
      )
      output[[3]]=tagList(radioButtons(
        inputId = "ImputationSelection",
        label = "Missing Value Imputation",
        choices = c("No Imputation" = "NA","Imputation" = "Imputation","Data without NAs"="CleanData"),
        selected = "NA"
      ),
      # Popovers for each radio button choice
     radioTooltip("ImputationSelection",choice = "NA",title = "No imputation on whole data"),
     radioTooltip("ImputationSelection",choice = "Imputation",title = "See documentation for details"),
     radioTooltip("ImputationSelection",choice = "CleanData",title = "No imputation and perform DE only on data without NAs")
      )
  }else{
    validate(
      need(ngroup !=2, "There are only two groups. More than three groups is required.")
    )
    
    output[[1]] = numericInput("count", "Number of Comparisons:", value = 3, min = 3)
    output$selectInputContainer <- renderUI({
      lapply(1:input$count, function(i) {
        selectInput(paste0("selectGroup", i), paste0("Select group", i, ":"), choices = vars, selected = vars[i])
      })
    })
    
  }
  output
})



#if submission button is clicked
observeEvent(input$DETestType,{
  progressSweetAlert(
    session = session,
    id = "DEProgress",
    title = "Read in data",
    display_pct = TRUE,
    value = 0
  )
  
  #based on inport load normalized data or original data if normalization not required
  if(is.data.frame(variables$corrected_data) &&nrow(variables$corrected_data)!=0){
    data <- variables$corrected_data
  }else if(is.data.frame(variables$normedData) &&nrow(variables$normedData)!=0)
    data <- variables$normedData
  else{
    data <- variables$CountData
  }
  
  
  updateProgressBar(
    session = session,
    id = "DEProgress",
    title = "Load Normalized Data",
    value = 10
  )
  
  #collect group information
  groupList <- variables$group
  group <- input$groups2
  idx <- as.numeric(grep(group, colnames(groupList)))
  if(colnames(groupList)[1] == "Sample"){
    colnames(groupList)[1] <- tolower(colnames(groupList)[1])
  }
  
  #get the new group
  variables$groupList <-
    lapply(unique(groupList[,idx]), function(x) {
      groupList[which(groupList[,idx] == x), ]$sample
    })
  names(variables$groupList) <- unique(groupList[,idx])
  data.cl <- rep(0, ncol(variables$CountData))
  
  for (i in 1:length(variables$groupList)) {
    data.cl[unlist(lapply(variables$groupList[[i]], convert2cl, df = variables$CountData))] <- names(variables$groupList[i])
  }
  
  data.cl <- data.cl[data.cl!=0]
  
  normed_tcc <- data
  metaData1 <- reactive(variables$group)
  data1 = reactive(list(data = as.data.frame(normed_tcc,),level = "protein"))
  
  updateProgressBar(
    session = session,
    id = "DEProgress",
    title = "Pre-processing Data",
    value = 20
  )
  
  #preprocess
  # input different treatments
  # Selection of a data subset (highly variable), This subset is for analyses
  df = data1()$data
  dfSample = metaData1()
  level = data1()$level
  metric = NULL
  pct = NULL
  res = preprocess(df, level, metric, pct)
  data2 = reactive(list(rawData = res$rawData, data = res$data, level = res$level, sampleInfo = dfSample))
  
  updateProgressBar(
    session = session,
    id = "DEProgress",
    title = "Retrieve group information",
    value = 30
  )

  # Differentially expressed peptides/proteins
  # Data processing
  df = data2()$data
  level = data2()$level
  dfSample = data2()$sampleInfo
  # Preparation of a statistical testing
  comparison = as.character()
  # factors = unique(dfSample[[input$groups2]])
  factors = c()
  if(input$compareSelection == "pair"){
    factors = c(input$group1select,input$group2select)
  }
  if(input$compareSelection == "multicomp"){
    for(i in 1:input$count){
      selectedValues <- input[[paste0("selectGroup", i)]]
      factors = c(factors,selectedValues)
    }
    #factors = c(input$group1select,input$group2select,input$group3select)
  }
  
  #Store factors as global variables
  variables$factors <- factors
  

  

  #factors = unique(factors[!is.na(factors)])
  nGroups = length(factors)
  
  for (g in 1:nGroups) {
    groupName = paste0("Group", g)
    #comparison[g] = paste(dfSample$ID[dfSample[[input$groups2]] == factors[g]], collapse = ",")    # Oldest version
    # comparison[g] = paste(dfSample[dfSample[[input$groups2]] == factors[g], 1], collapse = ",")
    cc = dfSample[dfSample[[input$groups2]] == factors[g], 1]
    cc <- unname(cc)
    cc <- unlist(cc)
    cc = cc[!is.na(cc)]
    comparison[g] = paste(cc, collapse = ",")
  }
  
  
  updateProgressBar(
    session = session,
    id = "DEProgress",
    title = "Imputation",
    value = 40
  )
  
  ############Imputation
  nGroups = length(comparison)
  groups = list()
  nSamples = 0
  samples = NULL
  for (g in 1:nGroups) {
    groups[[g]] = unlist(strsplit(comparison[g], ","))
    nSamples = nSamples + length(groups[[g]])
    samples[[g]] <- groups[[g]]
  }
  
  
  group1=df[,samples[[1]]]
  
  group2=df[,samples[[2]]]

  
  data_imputation <- cbind(group1,group2)
  
  
  #separate data into two groups
  count_non_na_per_group_base <- function(data, group1, group2) {
    count_group1 <- apply(group1, 1, function(x) sum(!is.na(x)))
    count_group2 <- apply(group2, 1, function(x) sum(!is.na(x)))
    result <- as.data.frame(cbind(count_group1,count_group2))
    return(result)
  }
  
  #perform imputation if selected
  # Count non-missing values for each group in each row
  data_counts <- count_non_na_per_group_base(data_imputation, group1, group2)
  
  
  
  impute_data <- function(data, group1, group2,counts) {
    
    for (i in 1:nrow(data)) {
      count1 <- counts$count_group1[i]
      count2 <- counts$count_group2[i]
      
      idx1 <- match(colnames(group1),colnames(data))
      
      idx2 <- match(colnames(group2),colnames(data))
      
      # Case 1: Only one group has more than 1 non-missing value
      if ((count1 > 1 && count2 < 1) || (count2 > 1 && count1 < 1)) {
        if (count1 > 1 && count2 < 1) {
          group_more <- group1
          idx_more <- idx1
          group_less <- group2
          idx_less <- idx2
        } else {
          group_more <- group2
          idx_more <- idx2
          group_less <- group1
          idx_less <- idx1
        }
        
        diff <- abs(count1-count2)
        diff <- min(diff,ncol(group_less),ncol(group_more))
        
        random_select_col <- sample(idx_less,diff)
        for (j in random_select_col) {
          data[i,j] <- min(data[j],na.rm = T)
        }


      }
      
      # Case 2: Both groups have less than 2 non-missing values
      if (count1 <= 1 && count2 <= 1) {
        data[i, ] <- NA  # Indicate to discard this row later
      }
      
      # Case 3: Both groups have more than 1 non-missing value
      # No action needed
    }
    
    # Discard rows with all NAs (case 2)
    data <- data[apply(data, 1, function(row) !all(is.na(row))), ]
    
    # Remove count columns before returning
    data <- data[, !(names(data) %in% c("count_group1", "count_group2"))]
    
    return(data)
  }
  
  
  df <- cbind(variables$CountData[,c(1,2)],df)
  #selection
  if (input$ImputationSelection == "Imputation") {
    data_imputed <- impute_data(df, group1, group2,data_counts)
  } else if (input$ImputationSelection == "CleanData") {
    data_imputed<-df[complete.cases(df[2:ncol(df)]),]
  } else {
    data_imputed<-df
  }
  
  description <- data_imputed[c(1,2)]
  data_imputed <- data_imputed[-c(1,2)]
  #perform limma
  statRes = reactive(statTest(data_imputed, level, comparison,dfSample,data.cl,factors))
  
  updateProgressBar(
    session = session,
    id = "DEProgress",
    title = "Calculate differentially expressed peptides/proteins",
    value = 50
  )
  


  # Show computation time notification
  #runtime <- round(normed_tcc$DEGES$execution.time[3], 2)
  
  updateProgressBar(
    session = session,
    id = "DEProgress",
    title = "Data Processing",
    value = 60
  )
  
  
  # Data processing
  shinyCatch({statres = statRes()})
  dfRaw = data2()$rawData
  exprs = statres$data
  sampleLabels = colnames(exprs)
  dfSample = data2()$sampleInfo
  # nGroups = length(unique(dfSample[[input$groups2]]))
  groupLabels = dfSample[[input$groups2]]
  groupLabels = groupLabels[!is.na(groupLabels)]
  #nGroups = length(unique(groupLabels))
  
  # Handle threshold inputs
  logFC = 0
  sigMetric = "p-value"
  sigCutoff = 1
  resLogFC = statres$res[, grep("Log2Fold", colnames(statres$res))]
  if (nGroups > 2) {
    absLogFC = apply(cbind(abs(apply(resLogFC, 1, min)), abs(apply(resLogFC, 1, max))), 1, max)
  } else {
    absLogFC = abs(resLogFC)
  }
  
  # Select DE peptides/proteins and organize a dataset for subsequent analyses
  rowInd = which(statres$res[[sigMetric]] < sigCutoff & absLogFC >= logFC)
  if (nGroups == 2) {
    exprs = cbind(description,exprs, `p-value` = statres$res$`p-value`, FDR = statres$res$FDR, Log2Fold = resLogFC)
  } else if (nGroups > 2) {
    exprs = cbind(description,exprs, `p-value` = statres$res$`p-value`, FDR = statres$res$FDR)
    exprs = cbind(exprs, resLogFC)
  }
  
  subData2 = reactive(list(rawData = dfRaw, sampleInfo = dfSample, data = exprs))
  
  
  res <- statRes()$res
  #res <- res[order(res$`p-value`),]
  result <- subData2()$data
  compareName <- NULL
  combMatrix = combn(seq(1, nGroups), 2)
  for (i in 1:ncol(combMatrix)) {
    compareName <- c(compareName, paste(factors[combMatrix[1, i]], "/", factors[combMatrix[2, i]], sep = ""))
  }
  ind <- grep("Log2Fold", colnames(result))
  for (k in 1:length(ind)) {
    colnames(result)[ind[k]] = paste("Log2Fold(", compareName[k], ")", sep = "")
  }

  
  variables$result <- result
  variables$res <- res
  
  output$downloadDE <- downloadHandler(
    filename = "DEresult.csv",
    content = function(file) {
      write.csv(result, file)
    }
  )
  
  # write.csv(variables$res, "test/res.csv")
  # write.csv(variables$result,"test/result.csv")
  
  # Render result table on the right top ----
  output$norm_resultTable <- DT::renderDataTable({
    if (nrow(variables$result) == 0) {
      DT::datatable(variables$result)
    } else {
       data = variables$result
       
       # colInd = grep('^sig', colnames(data))
       #data[, colInd] = formattable::comma(data[, colInd], digits = 3)
       data$`p-value` = formatC(data$`p-value`, digits = 4, format = "e")
       data$FDR = formatC(data$FDR, digits = 4,format = "e")
       log2Ind = grep("Log2Fold", colnames(data))
       colInd = c(3:(log2Ind[1]-3))
       data[colInd] = round(data[colInd], digits = 4)
       data[log2Ind] <- round(data[log2Ind],digits = 4)
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
  },server = T)
  
  closeSweetAlert(session = session)
  sendSweetAlert(
    session = session,
    title = "DONE",
    text = "Differential Experssion analysis done.",
    type = "success"
  )
  
  
  #store group information
  #get group information
  groupList <- as.data.frame(variables$group)
  group <- input$groups2
  idx <- as.numeric(grep(group, colnames(groupList)))
  
  factors <- variables$factors
  
  selected = groupList[groupList[,idx] %in% factors,]
  if(colnames(selected)[1] == "Sample"){
    colnames(selected)[1] <- tolower(colnames(selected)[1])
  }
  
  #get the new group
  variables$groupList2 <-
    lapply(unique(selected[,idx]), function(x) {
      selected[which(selected[,idx] == x), ]$sample
    })
  names(variables$groupList2) <- unique(selected[,idx])
  data.cl <- rep(0, ncol(variables$CountData))
  
  for (i in 1:length(variables$groupList2)) {
    data.cl[unlist(lapply(variables$groupList2[[i]], convert2cl, df = variables$CountData))] <- names(variables$groupList2[i])
  }
  data.cl <- data.cl[data.cl !=0]
  
  variables$data.cl_DE <- data.cl
  
  DERun$DERunValue <- input$DETestType
})


output$DEResultTable <- renderUI({
  if(DERun$DERunValue){
    tagList(
    fluidRow(column(
      12, 
      downloadButton("downloadDE", "Download Result Table"),
      DT::dataTableOutput('norm_resultTable') %>% withSpinner(),
      plotlyOutput('normboxplot')%>% withSpinner()
    )))} else {
      helpText("Click [Run Differential Analysis] to obtain Result Table.")
    }
})


observeEvent(input$norm_resultTable_rows_selected,{
  
  #input result data table after differential expression
  data <- variables$result
  #get only the values columns
  data = data %>% dplyr::select(-contains(c("Log2Fold","p-value","FDR")))
  colInd = c(3:ncol(data))
  # colInd = grep('^sig', colnames(data))
   data = data[, colInd]
  
  #get the row index when clicking the row
  rowInd = input$norm_resultTable_rows_selected
  
  #get group information
  groupList <- as.data.frame(variables$group)
  group <- input$groups2
  idx <- as.numeric(grep(group, colnames(groupList)))

  #factors <- variables$factors

  #selected = groupList[groupList[,idx] %in% factors,]

  if(colnames(groupList)[1] == "Sample"){
    colnames(groupList)[1] <- tolower(colnames(groupList)[1])
  }
  
  #get the new group
  variables$groupList <-
    lapply(unique(groupList[,idx]), function(x) {
      groupList[which(groupList[,idx] == x), ]$sample
    })
  
  names(variables$groupList) <- unique(groupList[,idx])
  
  data.cl <- rep(0, ncol(variables$CountData))
  
  for (i in 1:length(variables$groupList)) {
    data.cl[unlist(lapply(variables$groupList[[i]], convert2cl, df = variables$CountData))] <- names(variables$groupList[i])
  }
  
  data.cl <- data.cl[data.cl!=0]
  
  #Plot the bar plot when clicked on specific row
  output$normboxplot <- renderPlotly({
    expr = as.numeric(data[rowInd, ])
    df = data.frame(samples = colnames(data), intensity = round(expr,digits = 2))
    xOrder <-
      data.frame("name" = df$samples, "group" = data.cl)
    xOrderVector <- unique(xOrder[order(xOrder$group), ]$name)
    xform <- list(
      categoryorder = "array",
      categoryarray = xOrderVector,
      title = ""
    )
    
    plot_ly(
      data = df,
      x = ~samples,
      y = ~intensity,
      color = as.factor(data.cl),
      text = df$intensity,
      textposition = "outside",
      showlegend = FALSE,
      type = "bar",
      name = "DE_Result"
    )%>%
      layout(
        xaxis = xform,
        yaxis = list(title = "Intensity after differential analysis(log<sub>2)", range = list(0.9*min(df$intensity),max(df$intensity+0.5))),
        title = row.names(data[rowInd,])
      )%>%
      plotly::config(
        toImageButtonOptions = list(
          format = "svg",
          filename = row.names(data[rowInd,])
        )
      )
    
  })
  
  
})

#################Plots##############

########Moving SD################

output$DE_distributionUI <- renderUI({
  if(DERun$DERunValue && input$compareSelection == "pair"){
    fluidRow(column(
      3,
      numericInput(
        inputId = "DE_binsize",
        label = "Bin Size",
        value = 100
      ),
      do.call(actionBttn, c(
        list(
          inputId = "RunMovingSD",
          label = "Calculate moving SD of Log2FC",
          icon = icon("play")
        ),
        actionBttnParams
      ))
),
      column(
        9,
        plotlyOutput("DE_distributionPlot"),
        tags$hr(),
        downloadButton("downloadMSD", "Download MSD"),
        tags$br(),
        DT::dataTableOutput("MSDoutput")
      )
    )
  } else if(DERun$DERunValue && input$compareSelection == "multicomp"){
      helpText("Distribution plot is only available for pairwise comparison.")
    }
  else {
      helpText("Click [Run Differential Analysis] to obtain Result Table.")
    }
})



observeEvent(input$RunMovingSD,{
  #show notification of starting
  showNotification("Plotting Distribution...", type = "message")
  
  
  #retreive data
  data <- variables$result
  

  
  data = data %>% dplyr::select(-contains(c("p-value","FDR")))
  colInd = c(3:ncol(data))
  # colInd = grep('^sig', colnames(data))
  data = data[, colInd]
  
  n = ncol(data)
  
  
  #reorder data
  data$mean_intensity = rowMeans(data[-ncol(data)],na.rm = T)
  data <- data[order(data$mean_intensity,decreasing = F),]
  

  
  #get input binsize
  binsize <- input$DE_binsize
  
  # Calculate rolling standard deviation
  # rolling_sd <- rollapply(data$Log2Fold, 
  #                         width = binsize, 
  #                         by = 1, 
  #                         FUN = sd, na.rm = T,
  #                         align = 'right')
  
  # Calculate mean of the standard deviations of the next 100 proteins
  end <- (nrow(data))-binsize+1
  moving_sd <- sapply(1:end, function(i) {
    last <- i+binsize-1
    sd(data$Log2Fold[i:last],na.rm = T)
  })
  moving_sd_remain <- numeric(nrow(data) - end)
  moving_sd_remain[] <- moving_sd[length(moving_sd)]
  moving_sd <- c(moving_sd,moving_sd_remain)
  data$movingSD <- moving_sd
  
  data_row_names=rownames(data)
  
  data_reordered <- data[match(rownames(variables$result),data_row_names), ]
  data_MSD=cbind(data_reordered,variables$result[,c("p-value","FDR")])

  
  #generate distribution plot based on SD
  output$DE_distributionPlot <- renderPlotly({
    p <- plot_ly(
      data = data,
      x = ~mean_intensity,
      y = ~movingSD,
      type = "scatter",
      mode = "markers"
    )%>%
      layout(
        title = "Moving SD of log2 Fold Change",
        xaxis = list(title = "log2 protein intensity"),
        yaxis = list(title = "Moving SD",
                     range=c(0,1))
      )
    
    p
  })

  #output moving standard deviation table
  output$MSDoutput <- renderDataTable({
    if (nrow(variables$result) == 0) {
      DT::datatable(data_MSD)
    } else {
      data = data_MSD
      colInd <- ncol(data) - 5
      data[1:colInd]=round(data[1:colInd],digits = 4)
      log2Ind = grep("Log2Fold", colnames(data))
      data[, log2Ind] = round(data[, log2Ind], digits = 4)
      data$movingSD = round(data$movingSD,digits = 4)
      data$mean_intensity = round(data$mean_intensity,digits = 4)
      data$`p-value` = formatC(data$`p-value`, digits = 3, format = "e")
      data$FDR = formatC(data$FDR, digits = 3,format = "e")
      DT::datatable(
        data = data,
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
  })
  output$downloadMSD <- downloadHandler(
    filename = "MSD.csv",
    content = function(file) {
      write.csv(data_MSD, file)
    }
  )
  
  variables$DE_moving_sd <- data_MSD
})

################SD between group#############
output$DE_OldSDUI <- renderUI({
  if(DERun$DERunValue && input$compareSelection == "pair"){
    fluidRow(column(
      3,
      uiOutput("DE_distributionGroup"),
      uiOutput("DE_distributionGroup2")),
      column(
        9,
        plotOutput("DE_SDdistPlot")
      )
    )
  } else if(DERun$DERunValue && input$compareSelection == "multicomp"){
    helpText("Distribution plot is only available for pairwise comparison.")
  }
  else {
    helpText("Click [Run Differential Analysis] to obtain Result Table.")
  }
})


output$DE_distributionGroup <- renderUI({
  df = metaData1()
  if (is.null(colnames(df))) {
    vars = NULL
  } else {
    vars = colnames(df)[2: ncol(df)]    
  }
  
  selectInput("log2fcGroup", "Select Grouping variable to calculate log2FC cutoff",
              choices = vars, selected = vars[1])
  
})

output$DE_distributionGroup2 <- renderUI({
  group <- input$log2fcGroup
  output = tagList()
  df = metaData1()
  if (is.null(colnames(df)) || is.null(group)) {
    vars = NULL
  } else {
    idx = grep(group, colnames(df))
    vars = unique(df[,idx])
  }
  tagList(fluidRow(
    column(
      12,
      selectInput(
        inputId = "log2fcGroupselect",
        label = "Select variable to calculate SD of log2fc within this group",
        choices = vars,
        selected = vars[1]
      ),
      do.call(actionBttn, c(
        list(
          inputId = "RunDistribution",
          label = "Calculate SD of Log2FC",
          icon = icon("play")
        ),
        actionBttnParams
      ))
    )
  ))
  
  
})


observeEvent(input$RunDistribution,{
  #show notification of starting
  showNotification("Plotting Distribution...", type = "message")
  
  
  #retreive data
  data <- variables$result
  choice <- input$log2fcGroupselect
  group_select <- input$log2fcGroup
  
  data = data %>% dplyr::select(-contains(c("Log2Fold","p-value","FDR")))
  colInd = c(3:ncol(data))
  # colInd = grep('^sig', colnames(data))
  data = data[, colInd]
  
  #get group information
  group <- variables$groupList[choice]
  
  #select correlated samples
  selected_samples <- unlist(group)
  
  data <- data[selected_samples]
  n = ncol(data)
  
  
  #calculate log2fc for each sample
  log2fc = data
  for (i in 1:(ncol(data) - 1)) {
    for (j in (i + 1):ncol(data)) {
      log2fc[[paste0("log2FC(", names(data)[j], "/", names(data)[i],")")]] <- data[j] - data[[i]]
    }
  }
  
  log2fc = log2fc[-c(1:n)]
  
  #caluclate sd and mean for each sample log2fc
  sd <- apply(log2fc, 2, sd,na.rm = T)
  sd <- mean(sd)
  mean_logfc = apply(log2fc,2,mean,na.rm = T)
  mean_logfc <- mean(mean_logfc)
  
  
  #generate distribution plot based on SD and mean
  output$DE_SDdistPlot <- renderPlot({
    p <- ggplot()+
      stat_function(fun = dnorm,args = list(mean = mean_logfc, sd = sd))+
      theme_bw()+labs(x = "mean", y = "density")+xlim(-3,3)+
      annotate(geom = "text",
               label = paste0("SD = ",as.character(round(sd,5))),
               x = sd,
               y = 1,vjust = 1,hjust = 1.1)+
      annotate(geom = "text",
               label = paste0("2SD = ",as.character(round(2*sd,5))),
               x = 2*sd,
               y = 1,vjust = 1,hjust = -0.1)+
      geom_vline(xintercept=sd,lty=4,col="blue",lwd=0.4)+
      geom_vline(xintercept=2*sd,lty=4,col="blue",lwd=0.4)+
      ggtitle("Simulated distribution of SDs")
    
    p
  })
  
  
})



#Volcano plot
# Sample Distribution Density Plot UI ----
output$DE_volcanoPlot <- renderUI({
  if (DERun$DERunValue && input$compareSelection == "pair") {
    txt <- paste("(",input$group1select,"/",input$group2select,")")
    tagList(fluidRow(
      column(
        3,
        textInput("volplttitle", "Graphic Title", value = "Volcano Plot"),
        selectInput("metric2", 
                    label = "Select the measure of significance", 
                    choice = list("p-value" = "p-value", "FDR" = "FDR"), 
                    selected = 1),
        numericInput("pvalcutoff", 
                     label = "Significance level", 
                     min = 0, 
                     max = 1, 
                     step = 0.01, 
                     value = 0.05),
        numericInput("log2fcUp", 
                     label = "Upregulated Log2-fold cutoff", 
                     min = -2, 
                     max = 2, 
                     step = 0.1,
                     value = 0.1),
        numericInput("log2fcDown", 
                     label = "Downregulated Log2-fold cutoff", 
                     min = -2, 
                     max = 2, 
                     step = 0.1,
                     value = -0.1),
        do.call(actionBttn, c(
          list(
            inputId = "RunVolcano",
            label = "Run Volcano Plot",
            icon = icon("play")
          ),
          actionBttnParams
        )),
        tags$hr(),
        box(
          title = tagList(icon("option-vertical",lib = "glyphicon"), "Advanced Parameters"),
          solidHeader = F,
          width = NULL,
          collapsible = T,
          collapsed = T,
          textAreaInput(
            inputId = "labelProt",
            label = "Interested protein to label on plot",
            rows = 5,
            placeholder = "Input your protein accession number/Uniprot ID.\nOne protein each line"
          ),
          textInput(inputId = "volcxlab", 
                    "X-axis Label", 
                    value = paste("log<sub>2</sub>",txt)),
          sliderInput(
            "volcanoPointSize",
            "Point Size",
            min = 1,
            max = 8,
            value = 3,
            step = 0.2
          ),
          spectrumInput(
            inputId = "downColor",
            label = tagList("Down-regulated in G2", htmlOutput("downPreview")),
            choices = list(
              list(
                "#0000ff",
                "green",
                "black",
                "white",
                "blanchedalmond"
              ),
              as.list(brewer.pal(n = 9, name = "Blues")),
              as.list(brewer.pal(n = 9, name = "Greens")),
              as.list(brewer.pal(n = 11, name = "Spectral"))
            ),
            options = list(`toggle-palette-more-text` = "Show more")
          ),
          spectrumInput(
            inputId = "upColor",
            label = tagList("Up-regulated in G2", htmlOutput("upPreview")),
            choices = list(
              list(
                "red",
                "yellow",
                "orange",
                "white",
                "darkblue"
              ),
              as.list(brewer.pal(n = 9, name = "Reds")),
              as.list(brewer.pal(n = 9, name = "Oranges")),
              as.list(brewer.pal(n = 11, name = "Spectral"))
            ),
            options = list(`toggle-palette-more-text` = "Show more")
          )
        )
      ),
      column(
        9,
        plotlyOutput("volcanoPlot") %>% withSpinner(),
        tags$hr(),
        plotlyOutput("protBarPlotInVolcano") %>% withSpinner()
      )
    ))
  } else if(DERun$DERunValue && input$compareSelection == "multicomp"){
    helpText("Volcano plot is only available for pairwise comparison.")
  }
  else {
    helpText("No data for ploting. Please run differential expression analysis first.")
  }
})


#preview
observeEvent(
  {
    input$RunVolcano
  },
  {
    #Obtain the upper and lower log2fold change cutoff
    DownCut <- input$log2fcDown
    UpCut <- input$log2fcUp
    pValCut <- input$pvalcutoff
    sigmethod <- input$metric2
    
    #import the result data table
    res <- na.omit(variables$res)
    result <- variables$result
    result <- result[row.names(result) %in% row.names(res),]
    
    
    if(sigmethod == "p-value"){
      #get down part
      DownProt <- nrow(res[res$Log2Fold <= DownCut & res$`p-value` <= pValCut,])
      #get up part
      UpProt <- nrow(res[res$Log2Fold >= UpCut & res$`p-value` <= pValCut,])
    }else{
      DownProt <- nrow(res[res$Log2Fold <= DownCut & res$`FDR` <= pValCut,])
      UpProt <- nrow(res[res$Log2Fold >= UpCut & res$`FDR` <= pValCut,])
    }
    
    
    
    
    output$downPreview <- renderText({
      paste0(
        "<font color=\"",
        input$downColor,
        "\"><b>",
        DownProt,
        " proteins</b></font>"
      )
    })
    
    output$upPreview <- renderText({
      paste0(
        "<font color=\"",
        input$upColor,
        "\"><b>",
        UpProt,
        " proteins</b></font>"
      )
    })
    
    
    
    output$volcanoPlot <- renderPlotly({
      # Preparation of the statistical testing result for visualization
      # Check the number of groups for comparison, and move forward to the volcano plot
      if(dim(res)[2] > 3){
        sendSweetAlert(
          session = session,
          title = "ERROR",
          text = "Volcano Plot is unavailable for multiple comparison now.",
          type = "info"
        )
      }
      
      
      #input filter parameters
      sigMetric = input$metric2
      sigCutoff = input$pvalcutoff
      
      #remove the unnecessary column
      if (sigMetric == "p-value") {
        # res = res[, 2:3]
        res[, 3] = NULL
        ylab = "-log10(p-value)"
      } else if (sigMetric == "FDR") {
        # res = res[, c(2, 4)]
        res[, 2] = NULL
        ylab = "-log10(FDR)"
      }
      colnames(res) = c("logfc", "significance")
      
      res$color <- "None"
      
      if(DownProt!=0){
        res[which(res$logfc <= DownCut & res$significance <= sigCutoff),]$color <- "Down"
      }else{
        sendSweetAlert(
          session = session,
          type = "warning",
          title = "Warning",
          text = "No down regulated protein is filtered."
        )
      }
      if(UpProt!=0){
        res[which(res$logfc >= UpCut & res$significance <= sigCutoff),]$color <- "Up"
      }else{
        sendSweetAlert(
          session = session,
          type = "warning",
          title = "Warning",
          text = "No up regulated protein is filtered."
        )
      }

      
      level <- factor(res$color)
      levels(level) <- list(
        "Down" = 0,
        "None" = 1,
        "Up" = 2
      )
      
      #input volcano plot parameters
      xlab = input$volcxlab
      xmin = min(res$logfc)
      xmax = max(res$logfc)
      ymin = 0
      ymax = max(-log10(res$significance))
      
      #############Label User interested Proteins#############
      
      res$GN <- result$GN
      selectList <- row.names(res) %in% unlist(strsplit(x = input$labelProt, split = "[\r\n]"))
      top_prots <- res[selectList,]
      
      
      # # Find and label the top 10 proteins
      # res$GN <- result$GN
      # uprows <- res %>% filter(color %in% "Up")
      # top10_up <- head(uprows[order(uprows$logfc,decreasing = T),], 5)
      # downrows <- res %>% filter(color %in% "Down")
      # top10_down <- head(downrows[order(downrows$logfc),], 5)
      # top_prots <- rbind(top10_up,top10_down)
      
      # Add gene labels for all of the top genes we found
      # here we are creating an empty list, and filling it with entries for each row in the dataframe
      # each list entry is another list with named items that will be used by Plot.ly
      a <- list()
      for (i in seq_len(nrow(top_prots))) {
        m <- top_prots[i, ]
        a[[i]] <- list(
          x = m[["logfc"]],
          y = -log10(m[["significance"]]),
          text = m$GN,
          xref = "x",
          yref = "y",
          showarrow = TRUE,
          arrowhead = 4,
          arrowsize = 0.5,
          ax = 20,
          ay = -20
        )
      }

      #add annotation
      annotation <- row.names(res)
      
      # Volcano plot using plotly
      fig <- plot_ly(
        data = res,
        source = "volcano"
      )%>%add_trace(
        x = ~logfc, 
        y = ~-log10(significance),
        type = "scatter", 
        mode = "markers",
        color = ~level,
        colors = c(input$downColor, "black", input$upColor), 
        text = result$GN,
        marker = list(size = input$volcanoPointSize,
                      opacity = 0.7),
        hovertemplate = paste("ID: %{text}",
                              "<br>Log2fold: %{x:.4f}", 
                              "<br>sig. level: %{y}<extra></extra>"),
        key = ~annotation,
        showlegend = F
      )%>%layout(
        annotations = a
      )%>%plotly::config(
          toImageButtonOptions = list(
            format = "svg",
            filename = input$volplttitle
          )
        )
      
      fig = fig %>% add_segments(x = DownCut, xend = DownCut, y = ymin, yend = ymax, 
                                 line = list(dash = "dash", color = "black"),
                                 showlegend = FALSE)
      fig = fig %>% add_segments(x = UpCut, xend = UpCut, y = ymin, yend = ymax, 
                                 line = list(dash = "dash", color = "black"),
                                 showlegend = FALSE)
      fig = fig %>% add_segments(x = xmin - 0.5, xend = xmax + 0.5, y = -log10(sigCutoff), yend = -log10(sigCutoff), 
                                 line = list(dash = "dash", color = "black"),
                                 showlegend = FALSE)
      fig = fig %>% layout(yaxis = list(title = ylab,range = c(ymin, ymax + 0.5),
                                         linecolor = "rgba(0, 0, 0, 0)", linewidth = 0.5, ticks = "outside", mirror = TRUE),
                           xaxis = list(title = xlab,range = c(xmin - 0.5, xmax + 0.5),
                                        linecolor = "rgba(0, 0, 0, 0)", linewidth = 0.5, ticks = "outside", mirror = TRUE),
                           title = input$volplttitle,
                           showlegend = FALSE)
      variables$volcano <- fig
      fig
      
    })
    volrun$volrunval <- T
  }
  
)


#This function renders a bar plot when use mouse to click on each point
output$protBarPlotInVolcano <- renderPlotly({
  # Read in hover data
  eventdata <- event_data("plotly_click", source = "volcano")
  validate(need(
    !is.null(eventdata),
    "Click the point to show protein's signal level of interest."
  ))
  # Get point number
  prot_id <- eventdata$key
  # Get expression level (Original)
  expression <-
    variables$CountData[row.names(variables$CountData) == prot_id, ]
  # Get expression level (Normalized)
  expressionNor <-
    variables$result[row.names(variables$result) == prot_id, ]
  colInd = c(3:ncol(expressionNor))
  # colInd = grep('^sig', colnames(expressionNor))
  expressionNor <- expressionNor[,colInd]
  
  #get group list
  data.cl <- rep(0, ncol(variables$CountData))
  
  for (i in 1:length(variables$groupList2)) {
    data.cl[unlist(lapply(variables$groupList2[[i]], convert2cl, df = variables$CountData))] <- names(variables$groupList2[i])
  }
  data.cl <- data.cl[data.cl !=0]
  
  #order the expression data
  expression <- t(expressionNor[which(data.cl != 0)])
  expression <- round(expression,digits = 2)
  
  xOrder <-
    data.frame("name" = row.names(expression), "group" = data.cl)
  xOrderVector <- unique(xOrder[order(xOrder$group), ]$name)
  xform <- list(
    categoryorder = "array",
    categoryarray = xOrderVector,
    title = ""
  )
  
  plot_ly(
    x = ~ row.names(expression),
    y = ~ expression[, 1],
    color = as.factor(data.cl),
    text = expression[, 1],
    textposition = "outside",
    showlegend = FALSE,
    type = "bar",
    name = "Raw"
  )%>%
    layout(
      xaxis = xform,
      yaxis = list(title = "log<sub>2</sub> Intensity", range = list(0.9*min(expression[,1]),max(expression[,1])+0.5)),
      title = colnames(expression)
    ) %>%
    plotly::config(
      toImageButtonOptions = list(
        format = "svg",
        filename = colnames(expression)
      )
    )
  
  
})



#------------------HEATMAP-----------------
#Heatmap

output$DE_dendUI <- renderUI({
  if(DERun$DERunValue){
    tagList(fluidRow(
      column(
        3,
        uiOutput("heatmapParameter")
      ),
      column(
        9,
        if(HeatmapRun$heatmapRunValue){
          tagList(
            downloadButton("downloadHeatmap", "Download Heatmap Table"),
            plotlyOutput("heatmap", height = HeatmapRun$height)%>% withSpinner()
          )
        }
        else{
          helpText("No parameter to plot. Run Heatmap first.")
        }
      )
    ))
  }
  else{
    helpText("No data to plot. Run differential expression analysis first.")
  }
})

output$heatmapParameter <- renderUI({
  res <- variables$res
  tagList(
    radioGroupButtons(
      inputId = "heatmapProtSelectType",
      label = "Filter Differential Expression Result",
      choices = c(
        "List" = "By_list",
        "FDR" = "By_FDR",
        "Top" = "topPercent"
      ),
      justified = TRUE,
      status = "primary"
    ),
    uiOutput("heatmapSelectProt"),
    numericInput("heatmap_log2fc", 
                 label = "Log2-fold cutoff", 
                 min = 0, 
                 max = 10, 
                 step = 0.1,
                 value = 0),
    selectInput(
      "heatmapDist",
      "Distance Measure",
      choices = list(
        "Euclidean" = "euclidean",
        "Maximum" = "maximum",
        "Manhattan" = "manhattan",
        "Canberra" = "canberra",
        "Binary" = "binary",
        "Minkowski" = "minkowski"
      ),
      selected = "euclidean"
    ),
    selectInput(
      "heatmapCluster",
      "Agglomeration Method",
      choices = list(
        "ward.D" = "ward.D",
        "ward.D2" = "ward.D2",
        "Single" = "single",
        "Complete" = "complete",
        "UPGMA" = "average",
        "WPGMA" = "mcquitty",
        "WOGMC" = "median",
        "UPGMC" = "centroid"
      ),
      selected = "complete"
    ),
    selectInput(
      inputId = "dendCompute",
      label = "Dendrogram Selection",
      choices = c("Both" = "both",
                  "None" = "none",
                  "Row" = "row",
                  "Column" = "column"),
      selected = "both"
    ),
    selectInput(
      inputId = "colorSelectionMethod",
      label = "Color Selection Method",
      choices = c("Color map", "Two colors", "Three colors")
    ),
    uiOutput("heatmapColorSelectionPanel"),
    sliderInput(
      "heatmapColorNumber",
      "Select the number of colors to be in the palette",
      min = 1,
      max = 50,
      step = 1,
      value = 20
    ),
    tags$b("Color Preview"),
    plotOutput("colorPreview", height = "20px"),
    numericInput(
      inputId = "heatmapHeight",
      label = "Height of Heatmap",
      value = 1000,
      min = 500
    ),
    textInput(
      inputId = "heatmapTitle",
      label = "Heatmap Title",
      value = "Heatmap of top proteins"
    ),
    checkboxInput(
      inputId = "showrowlabel",
      label = "Show Row Label?",
      value = F),
    checkboxInput(
      inputId = "showdendrow",
      label = "Show Row Dendrogram?",
      value = T),
    
    do.call(actionBttn, c(
      list(
        inputId = "heatmapRun",
        label = "Run Heatmap",
        icon = icon("play")
      ),
      actionBttnParams
  )))
})


# Preview proteins count -----
observeEvent(input$heatmapFDR, {
  res <- variables$result
  prot_count <-
    nrow(res[res$FDR <= input$heatmapFDR, ])
  output$heatmapProteinPreview <- renderText({
    paste0(
      "Protein number: ",
      prot_count
    )
  })
})

output$heatmapSelectProt <- renderUI({
  res <- variables$res
  nprotein <- nrow(res)
  switch(
    input$heatmapProtSelectType,
    "By_list" = textAreaInput(
      "heatmapTextList",
      "Paste Protein List",
      rows = 5,
      placeholder = "Input uniprot ID (first column in the dataset), one protein id per line."
    ),
    "By_FDR" = tagList(
      tagList(
        sliderInput(
          "heatmapFDR",
          "FDR Cut-off",
          min = 0.01,
          max = 1,
          value = 0.01
        ),
        textOutput("heatmapProteinPreview")
      )
    ),
    "topPercent" = numericInput(
      inputId = "topProt",
      label = "Top Selected Proteins",
      value = 100,
      min = 0,
      max = nprotein,
      step = 10
    )
  )
    
})

observeEvent(input$colorSelectionMethod, {
  output$heatmapColorSelectionPanel <- renderUI({
    switch(input$colorSelectionMethod,
           "Color map" = {
             selectInput(
               "heatmapColor",
               "Choose Colormap",
               choices = list(
                 "PiYG", "PRGn", "BrBG", "PuOr", "OrRd",
                 "Oranges", "RdGy", "RdBu", "RdYlBu", "RdYlGn", 
                 "Spectral", "coolwarm"
               ),
               selected = "RdBu"
             )
           },
           "Two colors" = createColorSelectors(),
           "Three colors" = createColorSelectors(includeMiddle = TRUE)
    )
  })
})

createColorSelectors <- function(includeMiddle = FALSE) {
  tagList(
    spectrumInput(
      inputId = "heatmapTwoColorLow",
      label = "Low",
      choices = list(
        list(
          "blue", "black", "gray30", "white", "blanchedalmond",
          "#2066a8", "#3594cc", "forestgreen"
        ),
        as.list(brewer.pal(n = 9, name = "Blues")),
        as.list(brewer.pal(n = 9, name = "Greens")),
        as.list(brewer.pal(n = 11, name = "Spectral")),
        as.list(brewer.pal(n = 8, name = "Dark2"))
      ),
      options = list(`toggle-palette-more-text` = "Show more")
    ),
    if (includeMiddle) {
      spectrumInput(
        inputId = "heatmapTwoColorMiddle",
        label = "Middle",
        choices = list(
          list(
            "white", "gray50", "blanchedalmond", "steelblue", 
            "#0000FF", "forestgreen"
          ),
          as.list(brewer.pal(n = 9, name = "Blues")),
          as.list(brewer.pal(n = 9, name = "Greens")),
          as.list(brewer.pal(n = 11, name = "Spectral")),
          as.list(brewer.pal(n = 8, name = "Dark2"))
        ),
        options = list(`toggle-palette-more-text` = "Show more")
      )
    },
    spectrumInput(
      inputId = "heatmapTwoColorHigh",
      label = "High",
      choices = list(
        list(
          "red", "brown", "orange", "white", "steelblue", 
          "forestgreen"
        ),
        as.list(brewer.pal(n = 9, name = "Blues")),
        as.list(brewer.pal(n = 9, name = "Greens")),
        as.list(brewer.pal(n = 11, name = "Spectral")),
        as.list(brewer.pal(n = 8, name = "Dark2"))
      ),
      options = list(`toggle-palette-more-text` = "Show more")
    )
  )
}

colorPanel <- reactive({
  if (input$colorSelectionMethod == "Color map" && length(input$heatmapColor) > 0) {
    colorPal <- switch(
      input$heatmapColor,
      "PiYG" = rev(PiYG(input$heatmapColorNumber)),
      "PRGn" = PRGn(input$heatmapColorNumber),
      "BrBG" = BrBG(input$heatmapColorNumber),
      "PuOr" = PuOr(input$heatmapColorNumber),
      "OrRd" = OrRd(input$heatmapColorNumber),
      "Oranges" = Oranges(input$heatmapColorNumber),
      "RdGy" = rev(RdGy(input$heatmapColorNumber)),
      "RdBu" = rev(RdBu(input$heatmapColorNumber)),
      "RdYlBu" = rev(RdYlBu(input$heatmapColorNumber)),
      "RdYlGn" = rev(RdYlGn(input$heatmapColorNumber)),
      "Spectral" = Spectral(input$heatmapColorNumber),
      "coolwarm" = cool_warm(input$heatmapColorNumber)
    )
  } else if (input$colorSelectionMethod == "Two colors" && length(input$heatmapTwoColorLow) > 0) {
    colorPal <- colorRampPalette(c(input$heatmapTwoColorLow, input$heatmapTwoColorHigh))(input$heatmapColorNumber)
  } else if (input$colorSelectionMethod == "Three colors" && length(input$heatmapTwoColorLow) > 0) {
    colorPal <- colorRampPalette(c(input$heatmapTwoColorLow, input$heatmapTwoColorMiddle, input$heatmapTwoColorHigh))(input$heatmapColorNumber)
  } else {
    colorPal <- c("white")
  }
  
  colorPal
})

output$colorPreview <- renderPlot({
  colorPal <- colorPanel()
  op <- par(mar = c(0.5, 0, 0, 0))
  plot(
    c(0, length(colorPal)),
    c(0, 1),
    type = "n",
    xlab = "",
    ylab = "",
    ann = FALSE,
    bty = "n",
    xaxt = "n",
    yaxt = "n"
  )
  rect(0:(length(colorPal) - 1), 0, 1:(length(colorPal)), 1, col = colorPal, lwd = 0)
  par(op)
})



observeEvent(input$heatmapRun,{
  progressSweetAlert(
    session = session,
    id = "heatmapProgress",
    title = "Work in progress",
    display_pct = TRUE,
    value = 0
  )
  
  updateProgressBar(
    session = session,
    id = "heatmapProgress",
    title = "Processing data",
    value = 10
  )
  
  data <- variables$result
  res <- variables$res
  dfSample = variables$sampleInfo
  colorPal <- colorPanel()
  
  updateProgressBar(
    session = session,
    id = "heatmapProgress",
    title = "Processing selected lists",
    value = 30
  )
  
  #get the selected number of proteins for analysis
  if(input$heatmapProtSelectType == "By_list"){
    selectList <- row.names(data) %in% unlist(strsplit(x = input$heatmapTextList, split = "[\r\n]"))
    sampleData <- data[selectList,]
    
    updateTextAreaInput(session = session, inputId = "heatmapTextList",value = input$heatmapTextList)
  }
  else if(input$heatmapProtSelectType == "By_FDR"){
    selectList <-
      row.names(data) %in% row.names(res[res$FDR <= input$heatmapFDR, ])
    sampleData <- data[selectList,]
    
    updateSliderInput(session = session,inputId = "heatmapFDR",value = input$heatmapFDR)
  }
  else if(input$heatmapProtSelectType == "topPercent"){
    topProt <- input$topProt
    tempdata <- data[order(data$`p-value`), ]
    sampleData <- tempdata[1:topProt,]
    
    updateNumericInput(session = session,inputId = "topProt",value = input$topProt)
  }
  else{
    sendSweetAlert(
      session = session,
      title = "ERROR",
      text = "Select a method to plot heatmap first!",
      type = "info"
    )
    return()
  }
  updateProgressBar(
    session = session,
    id = "heatmapProgress",
    title = "Updating parameters",
    value = 50
  )
  
  #filter with Log2Fold threshold
  log2f <- input$heatmap_log2fc
  log2filterUp <- res %>%
    dplyr::select(contains("Log2Fold")) %>%
    filter(if_any(contains("Log2Fold"), ~ . >= log2f))
  
  selectList <- row.names(sampleData) %in% row.names(log2filterUp)
  df1 <- sampleData[selectList, ]
  
  log2filterDown <- res %>%
    dplyr::select(contains("Log2Fold")) %>%
    filter(if_any(contains("Log2Fold"), ~ . <= -log2f))
  selectList <- row.names(sampleData) %in% row.names(log2filterDown)
  df2 <- sampleData[selectList,]
  
  
  #merge datasets
  sampleData <- merge(df1,df2,by = "row.names", all = TRUE, suffixes = c(".Up", ".Down"))
  rownames(sampleData) <- sampleData$Row.names
  sampleData <- sampleData[ , -1] # Remove the Row.names column
  sampleData <- sampleData[1:ncol(data)]
  colnames(sampleData) <- colnames(data)
  
  
  #update input values
  updateNumericInput(session = session,inputId = "heatmap_log2fc",value = input$heatmap_log2fc)
  updateNumericInput(session = session, inputId = "heatmapHeight", value = input$heatmapHeight)
  updateTextInput(session = session, inputId = "heatmapTitle", value = input$heatmapTitle)
  
  #remove unnecessary columns
  sampleData = sampleData %>% dplyr::select(-contains(c("Log2Fold","p-value","FDR")))
  sampleData = sampleData[3:ncol(sampleData)]

  #change to heatmap format
  mat = as.matrix(sampleData)
  
  #check NAs
  giveNAs = which(is.na(as.matrix(dist(mat))),arr.ind=TRUE)
  tab = sort(table(c(giveNAs)),decreasing=TRUE)
  checkNA = sapply(1:length(tab),function(i){
    sum(is.na(as.matrix(dist(mat[-as.numeric(names(tab[1:i])),]))))
  })
  rmv = names(tab)[1:min(which(checkNA==0))]
  
  mat <- mat[complete.cases(mat),]
  mat = mat[-as.numeric(rmv),]
  
  mat = t(scale(t(mat), center = T, scale = F)) # Only mean-centering
  
  output$downloadHeatmap <- downloadHandler(
    filename = "heatmap.csv",
    content = function(file) {
      write.csv(mat, file)
    }
  )
  
  
  #plot heatmap
  output$heatmap <- renderPlotly({
    isolate(HeatmapRun$height <- input$heatmapHeight)
    p <- heatmaply(
      mat,
      #k_col = length(variables$factors),
      colors = colorPal,
      na.value = "grey50",
      na.rm = TRUE,
      dist_method = input$heatmapDist,
      hclust_method = input$heatmapCluster,
      xlab = "Sample",
      ylab = "Protein ID",
      dendrogram = input$dendCompute,
      showticklabels=c(TRUE, input$showrowlabel),
      show_dendrogram = c(input$showdendrow, T),
      main = input$heatmapTitle,
      key.title = "z score"
      )%>%
      plotly::config(
        toImageButtonOptions = list(
          format = "svg",
          filename = input$heatmapTitle
        ))
    
    variables$de_heatmap <- p
    p
  })
  
  
  updateProgressBar(
    session = session,
    id = "heatmapProgress",
    title = "Extracting data",
    value = 90
  )
  updateProgressBar(
    session = session,
    id = "heatmapProgress",
    title = "All done",
    value = 100
  )
  HeatmapRun$heatmapRunValue <- input$heatmapRun
  closeSweetAlert(session = session)
  sendSweetAlert(
    session = session,
    title = "Completed!",
    type = "success"
  )
  
  
})
