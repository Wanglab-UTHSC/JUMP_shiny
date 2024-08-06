# server-data-import.R

options(shiny.maxRequestSize = 60*1024^2)

dataImportCheck <- reactiveValues(importRunValue = FALSE)

output$dataSourceSelect <- renderUI({
  selectInput(
    "SampleDatabase",
    "Select Sample Data",
    choices = c(
      "JUMP sample" = "sample_data/combined_norm_uni_test_example.txt"
    )
  )
})

# If Sample Data button has been clicked, load sample data. ----

observeEvent(input$CountDataSample, {
  dataImportCheck$importRunValue <- FALSE

  if (input$SampleDatabase == "sample_data/combined_norm_uni_test_example.txt") {
    # data(hypoData)
    variables$CountData <- data.frame(fread(input$SampleDatabase), row.names = 1)
  } else {
    variables$CountData <- data.frame(variables$simulationData$count)
  }
  
  sendSweetAlert(
    session = session,
    title = "DONE",
    text = "Count data were successfully loaded.",
    type = "success"
  )
  
})


# If Upload Data button has been clicked, load the data via upload ----

dupValidate <- function(input){
  rownum <- nrow(input[duplicated(input$pid),])
  if (rownum > 0) {
    sendSweetAlert(
      session = session,
      title = "Duplicated accession numbers found!",
      text = "Unique value are required. Duplicated accession number rows are automatically dropped.",
      type = "warning"
    )
    return()
  } else {
    NULL
  }
}

observeEvent(input$uploadCountData, {
  showNotification("Start uploading file...", type = "message")
  
  shinyCatch(
    {
      #read in uploaded data
      rawdata <-
        data.frame(fread(input$uploadCountData$datapath))
      
      #omit NA values
      nNA <- sum(is.na(rawdata))
      variables$nNA <- nNA
      #rawdata<- na.omit(rawdata)
      
      #check the format
      if(input$dataType == "jumpq"){
        if(ncol(rawdata) <=24) stop("Input data error! Please check your selected data type and format.")
        colnum <- c(2:4,24:ncol(rawdata))
        rawdata <- rawdata[,colnum]
        colnames(rawdata)[1:3] <- c("pid","Description","GN")
        
      }
      else if(input$dataType  == "jump_batch"){
        batchidx <- grep("^batch",colnames(rawdata))
        
        #check if "Batch" is present in sample data name
        if(length(batchidx) == 0){
          sendSweetAlert(
            session = session,
            title = "Batch Data Error",
            text = "Data should have `batch` in sample names or change selected data type.",
            type = "error"
          )
          return()
        }
        colnum <- 2:4
        for(i in 1:length(batchidx)){
          j = i+1
          if(j <= length(batchidx)){
            first <- batchidx[i]+5
            last <- batchidx[j]-1
            colnum <- c(colnum, first:last)
          }
          
        }
        first <- batchidx[-1]+5
        last <- ncol(rawdata)
        colnum <- c(colnum, first:last)
        rawdata <- rawdata[,colnum]
        colnames(rawdata)[1:3] <- c("pid","Description","GN")
      }
      else{
        colnames(rawdata)[1:3] <- c("pid","GN","Description")
      }
    }
  )   
      
      
      #Sort raw data based on summation of signals
      # rawdata$sum <- rowSums(rawdata[,c(4:ncol(rawdata))])
      # rawdata <- rawdata[order(rawdata$sum,decreasing = T),]
      # 
      # 
      # rawdata<- rawdata[,-ncol(rawdata)]
      new <- rawdata
      
      
      tryCatch({
        #make the first column as row name
        rn <- rawdata$'pid'
        
        #check if the accession number is already obtained
        if(sum(grepl('\\|', rn))){
          rawdata2 <- rawdata %>% separate("pid", into = c("id1", "pid","id2"), sep = "\\|")
          pid = rawdata2$pid
          new <- data.frame(pid, rawdata)
          new <- new[-2]
        }
      },
      error = function(e){
        sendSweetAlert(
          session = session,
          title = "Accession number or protein id error. Check your first column.",
          text = as.character(message(e)),
          type = "error"
        )
        return()
      })
      
      
      #Get unique row names
      validate(dupValidate(new))
      new <- new[!duplicated(new$pid),]
      variables$dup <- nrow(new[duplicated(new$pid),])
      row.names(new) <- new$pid
      new <- new[,-1]
      
      
      #save the raw data as a global variable
      variables$CountData <- new

      dataImportCheck$importRunValue <- FALSE

      showNotification("Received uploaded file.", type = "message")

})

datasetInput <- reactive({
  variables$CountData
})


# Render a table of raw count data, adding color ----

output$table <- DT::renderDataTable({
  df <- datasetInput()
  temp <- df[,3:ncol(variables$CountData)]
  temp <- round(temp,digits = 0)
  df <- cbind(variables$CountData[,c(1,2)],temp)
  
  # brks <-
  #   quantile(df %>% select_if(is.numeric),
  #            probs = seq(.05, .95, .05),
  #            na.rm = TRUE
  #   )
  
  DT::datatable(
    df,
    colnames = c("Accesion Number" = 1),
    extensions = c("Scroller","RowReorder"),
    option = list(
      rowReorder = TRUE,
      deferRender = TRUE,
      autoWidth = TRUE,
      scrollY = 400,
      scroller = TRUE,
      scrollX = TRUE,
      pageLength = 5,
      searchHighlight = TRUE,
      orderClasses = TRUE
    )
  ) 
},server = T)


# Render DataTable of row data count ----

output$emptyTable <- renderUI({
  if (nrow(datasetInput()) == 0) {
    tags$p("No data to show. Click", tags$code("Sample"), "or", tags$code("Upload"), "your own dataset.")
  } else {
    DT::dataTableOutput("table")
  }
})

observeEvent(input$confirmedGroupList, {
  if (nrow(datasetInput()) == 0) {
    sendSweetAlert(
      session = session,
      title = "ERROR",
      text = "Please input raw data table!",
      type = "error"
    )
    return()
  }
  if (is.null(input$uploadGroup)) {
    sendSweetAlert(
      session = session,
      title = "ERROR",
      text = "Please input group information!",
      type = "error"
    )
    return()
  }
  
  tryCatch(
    {
      progressSweetAlert(
        session = session,
        id = "dataImportProgress",
        title = "Processing group info",
        display_pct = TRUE,
        value = 0
      )
      showNotification("Start uploading file...", type = "message")
      
      
      group <- data.frame(fread(input$uploadGroup$datapath, header = T, blank.lines.skip = T))
      group_no_header <- data.frame(fread(input$uploadGroup$datapath,header = FALSE))
      
      showNotification("Received uploaded file.", type = "message")
      
      #store group information
      variables$group <- group
      
      
      # variables$group_no_header <- group_no_header
      group_no_header <- group_no_header[-1,]
      groupName <- colnames(group)
      
      
      
      variables$groupList <-
        lapply(unique(group_no_header$V2), function(x) {
          group_no_header[group_no_header$V2 == x, ]$V1
        })
      names(variables$groupList) <- unique(group_no_header$V2)
      
      data.cl <- rep(0, ncol(variables$CountData))
      
      for (i in 1:length(variables$groupList)) {
        data.cl[unlist(lapply(variables$groupList[[i]], convert2cl, df = variables$CountData))] <- names(variables$groupList[i])
      }
      
      # Storage convert group list to local
      variables$groupListConvert <- data.cl
      
      #store numeric data
      variables$count.data <- variables$CountData[data.cl != 0]
      variables$count.data <- as.matrix(variables$count.data)
      
      #check whether sample info mathes raw data
      colnum <- ncol(variables$count.data)
      sample_nrow <- nrow(group)
      
      if(colnum != sample_nrow){
        sendSweetAlert(
          session = session,
          title = "Warning",
          text = "Sample number should match group info.",
          type = "warning"
        )
        return()
      }
      
      #store converted group info
      group2 <- data.cl[data.cl != 0]
      group2 <- data.frame(group = group2)
      group2 <- data.frame(lapply(group2, as.factor))
      rownames(group2) <- colnames(variables$count.data)
      variables$group_import <- group2
      
      updateProgressBar(
        session = session,
        id = "dataImportProgress",
        title = "Summarizing data",
        value = 90
      )
      
      output$groupSelection = renderUI({
        
        meatData1 <- reactive(variables$group)
        df = metaData1()
        if (is.null(colnames(df))) {
          vars = NULL
        } else {
          vars = colnames(df)[2: ncol(df)]    
        }

        selectInput("groups1", "Grouping variable",
                    choices = vars, selected = vars[1])
      })
      
      closeSweetAlert(session = session)
      sendSweetAlert(
        session = session,
        title = "DONE",
        text = "Group labels were successfully assigned.",
        type = "success"
      )
      
      dataImportCheck$importRunValue <- input$confirmedGroupList
      
    },
    error = function(e) {
      sendSweetAlert(
        session = session,
        title = "ERROR",
        text = "Check your group information format!",
        type = "error"
      )
      return()
    },
    warning = function(w) {
      sendSweetAlert(
        session = session,
        title = "Group Error!",
        text = "Check your group information format!",
        type = "error"
      )
      return()
    }
  )
})


# 
observeEvent(input$groups1,{
  #import the whole group selection
  groupList <- as.data.frame(variables$group)
  group <- input$groups1

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
  
  # Storage convert group list to local
  variables$groupListConvert <- data.cl
  

  variables$count.data <- variables$CountData[data.cl != 0]
  variables$count.data <- as.matrix(variables$count.data)
  
  #store converted group info
  group2 <- data.cl[data.cl != 0]
  group2 <- data.frame(group = group2)
  group2 <- data.frame(lapply(group2, as.factor))
  rownames(group2) <- colnames(variables$count.data)
  variables$group_import <- group2


})

output$importDataSummary <- renderUI({
  dt <- datasetInput()
  dup <- variables$dup
  
  rowCount <- nrow(dt)
  groupCount <- length(variables$groupList)
  groupText <- sapply(variables$groupList, length)
  if (length(groupText) > 0) {
    gText <- paste0(names(groupText), ": ", groupText, collapse = "\n")
  } else {
    gText <- NULL
  }
  variables$groupText <- gText
  
  data <- variables$CountData
  data.cl <- variables$groupListConvert
  cName <- unlist(variables$groupList)
  nNA <- variables$nNA
  
  tagList(
    tipify(
      tags$p(tags$b("N", tags$sub("Protein")), ":", rowCount),
      title = "Number of Proteins",
      placement = "left"
    ),
    tipify(
      tags$p(tags$b("N", tags$sub("group")), ": ", groupCount),
      title = "Number of Groups",
      placement = "left",
    ),
    tipify(
      tags$p(tags$b("N", tags$sub("samples in group")), ": ", gText),
      title = "Number of samples in each group",
      placement = "left"
    ),
    tipify(
      tags$p(tags$b("N", tags$sub("replicates")), ": ", dup),
      title = "Number of replicated gene names",
      placement = "left"
    ),
    tipify(
      tags$p(tags$b("N", tags$sub("NAs")), ": ", nNA),
      title = "Number of NAs omitted",
      placement = "left"
    )
  )
})

v <- reactiveValues(importActionValue = FALSE)

# This function render a boxplot of sample distribution ----

output$sampleDistributionBox <- renderPlotly({
  if (length(variables$count.data) > 0) {

    data <- variables$count.data
  
    
    log2data <- as.data.frame(log2(data))
    
    
    # Filter
    log2data <- log2data %>% filter_all(all_vars(. > input$sampleDistributionFilterLow | is.na(.)))
    
    data_stack <- data.frame(stack(log2data))
    
    # Add Group info
    group <-
      data.frame(
        "ind" = rownames(variables$group_import),
        "group" = variables$group_import$group
      )
    
    data <- left_join(data_stack, group, by = "ind")
    data <- arrange(data, group)
    
    p <- plot_ly(
      data = data,
      x = ~ind,
      y = ~values,
      type = "box",
      split = ~group,
      color = ~group
    ) %>%
      layout(
        title = input$sampleDistributionTitle,
        xaxis = list(title = input$sampleDistributionXlab, categoryarray = "array", categoryarray = ~col),
        yaxis = list(title = input$sampleDistributionYlab)
      ) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "svg",
          filename = input$sampleDistributionTitle
        )
      )
    variables$sampleDistributionBar <- p
    p
  } else {
    return()
  }
})

# Render a density plot of sample distribution ----
output$sampleDistributionDensity <- renderPlotly({
  if (length(variables$count.data) > 0) {

    dt <- as.data.frame(variables$count.data)
    if (input$densityFilter != "Do not filter") {
      count <-
        dt %>% filter_all(all_vars(. > input$densityFilter| is.na(.)))
    } else {
      count <- dt
    }
    data <- log2(count)
    
    group <- variables$group_import

    densityTable <- apply(data, 2, function(x) {
      density(x)
    })
    p <- plot_ly(type = "scatter", mode = "lines")
    for (i in 1:length(densityTable)) {
      p <- add_trace(
        p,
        x = densityTable[[i]][[1]],
        y = densityTable[[i]][[2]],
        color = group[rownames(group) == names(densityTable[i]), ],
        name = names(densityTable[i])
      )
    }
    
    pp <- p %>%
      layout(
        title = input$sampleDistributionDenstityTitle,
        xaxis = list(title = input$sampleDistributionDensityXlab),
        yaxis = list(title = input$sampleDistributionDensityYlab)
      ) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "svg",
          filename = input$sampleDistributionDenstityTitle
        )
      )
    variables$sampleDistributionDensity <- pp
    pp
  } else {
    return()
  }
})

# Sample Distribution Density Plot UI ----
output$sampleDistributionDensityPanel <- renderUI({
  if (dataImportCheck$importRunValue) {
    tagList(fluidRow(
      column(
        3,
        popify(
          helpText("Filter proteins with a total read count smaller than thresholds."),
          title = "Reference",
          content = 'Sultan, Marc, et al. <a href="http://science.sciencemag.org/content/321/5891/956">"A global view of gene activity and alternative splicing by deep sequencing of the human transcriptome."</a> <i>Science</i> 321.5891 (2008): 956-960.',
          placement = "left"
        ),
        sliderTextInput(
          inputId = "densityFilter",
          label = "Filter proteins threshold",
          choices = c("Do not filter", c(0:30))
        ),
        textInput(
          inputId = "sampleDistributionDenstityTitle",
          label = "Title",
          value = "Raw Intensity",
          placeholder = "Original Raw Intensity"
        ),
        textInput(
          inputId = "sampleDistributionDensityXlab",
          label = "X label",
          value = "log<sub>2</sub>(Intensity)",
          placeholder = "log<sub>2</sub>(Intensity)"
        ),
        textInput(
          inputId = "sampleDistributionDensityYlab",
          label = "Y label",
          value = "Density",
          placeholder = "Density"
        )
      ),
      column(
        9,
        plotlyOutput("sampleDistributionDensity") %>% withSpinner()
      )
    ))
  } else {
    helpText("No data for ploting. Please import dataset and assign group information first.")
  }
})

# Sample Distribution Boxplot UI ----
output$sampleDistributionBoxPanel <- renderUI({
  if (dataImportCheck$importRunValue) {
    tagList(fluidRow(
      column(
        3,
        sliderInput(
          inputId = "sampleDistributionFilterLow",
          label = "Filter low intensity proteins",
          min = 0,
          max = 30,
          value = 0
        ),
        textInput(
          inputId = "sampleDistributionTitle",
          label = "Title",
          value = "Raw Intensity",
          placeholder = "Raw Intensity"
        ),
        textInput(
          inputId = "sampleDistributionXlab",
          label = "X label",
          value = "Sample",
          placeholder = "Sample"
        ),
        textInput(
          inputId = "sampleDistributionYlab",
          label = "Y label",
          value = "log<sub>2</sub>(Intensity)",
          placeholder = "log<sub>2</sub>(Intensity)"
        )
      ),
      column(
        9,
        plotlyOutput("sampleDistributionBox") %>% withSpinner()
      )
    ))
  } else {
    helpText("No data for ploting. Please import dataset and assign group information first.")
  }
})

# Filtering low count under different low number, barplot ----
output$lowCountFilterByCutoff <- renderPlotly({
  if (length(variables$count.data) > 0) {

    data <- variables$count.data
    originalCount <- nrow(data)
    
    lowCount <-
      sapply(0:input$lowCountSlide, function(x) {
        sum(rowSums(data) > x)
      })
    
    lowCountdt <- data.frame(
      "Cutoff" = 0:input$lowCountSlide,
      "Filtered" = originalCount - lowCount,
      "Remain" = lowCount
    )
    
    plot_ly(
      lowCountdt,
      name = "Remain",
      x = ~Cutoff,
      y = ~Remain,
      text = ~Remain,
      textposition = "outside",
      hoverinfo = "text+name",
      hovertext = ~ paste0(
        "Cut off: ",
        Cutoff,
        "<br>Filtered number: ",
        Filtered,
        "<br>Remain number: ",
        Remain,
        "(",
        round(Remain / nrow(data) * 100, 2),
        "%)"
      ),
      type = "bar"
    ) %>%
      add_trace(
        name = "Filtered",
        y = ~Filtered,
        text = ~Filtered,
        type = "bar",
        textposition = "inside"
      ) %>%
      layout(
        title = "Filtering Threshold for Low Intensity",
        barmode = "stack",
        xaxis = list(title = "Filtering Low Intensity Cut off"),
        yaxis = list(title = "Protein number")
      ) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "svg",
          filename = "Filtering_Threshold_for_Low_Count_Proteins"
        )
      )
  } else {
    return()
  }
})

# Render Filtering Cutoff UI ----
output$lowCountFilterByCutoffUI <- renderUI({
  if (dataImportCheck$importRunValue) {
    tagList(fluidRow(
      column(
        3,
        popify(
          helpText("Filter proteins with a total read count smaller than thresholds."),
          title = "Reference",
          content = 'Sultan, Marc, et al. <a href="http://science.sciencemag.org/content/321/5891/956">"A global view of protein activity and alternative splicing by deep sequencing of the human transcriptome."</a> <i>Science</i> 321.5891 (2008): 956-960.',
          placement = "left"
        ),
        sliderInput(
          inputId = "lowCountSlide",
          label = "Max threshold",
          min = 3,
          max = 50,
          value = 15,
          step = 1
        )
      ),
      column(
        9,
        plotlyOutput("lowCountFilterByCutoff") %>% withSpinner()
      )
    ))
  } else {
    helpText("No data for ploting. Please import dataset and assign group information first.")
  }
})

# MDS Plot ----

output$mdsPlotObject <- renderPlotly({
  if (length(variables$count.data) > 0) {

    dt <- variables$count.data
    if (input$mds == "Nonmetric MDS") {
      mds <-
        data.frame(isoMDS(dist(
          1 - cor(dt, method = input$mdsMethod),
          method = input$mdsDistMethod
        )))
    } else {
      mds <-
        data.frame(cmdscale(dist(
          1 - cor(dt, method = input$mdsMethod),
          method = input$mdsDistMethod
        )))
    }
    mds$name <- rownames(mds)
    mdsG <- variables$group_import
    mdsG$name <- rownames(mdsG)
    mdsJ <- left_join(mds, mdsG, by = "name")
    p <- plot_ly(
      data = mdsJ,
      x = mdsJ[, 1],
      y = mdsJ[, 2],
      type = "scatter",
      mode = "text",
      text = ~name,
      color = ~group
    ) %>%
      layout(title = paste0(input$mds, " Plot")) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "svg",
          filename = paste0(input$mds, "_Plot")
        )
      )
    
    variables$mdsPlotplot <- p
    variables$mdsPlot[["params"]] <- list(
      "mds" = input$mds,
      "mdsMethod" = input$mdsMethod,
      "mdsDistMethod" = input$mdsDistMethod
    )
    p
  } else {
    return()
  }
})

# MDS information
output$MDShelpText <- renderUI({
  helpText(
    "Use all proteins' raw intensity to calculate",
    input$mdsMethod,
    "correlation coefficient (rho) to create a matrix of (1 - rho). Calculate the ",
    input$mdsDistMethod,
    " distances between samples and plot the result to a two-dimension MDS plot."
  )
})

# Render MDS plot ----
output$mdsUI <- renderUI({
  if (dataImportCheck$importRunValue) {
    tagList(fluidRow(
      column(
        3,
        uiOutput("MDShelpText"),
        selectInput(
          inputId = "mdsMethod",
          label = "Correlation Coefficient",
          choices = c(
            "Spearman" = "spearman",
            "Pearson" = "pearson",
            "Kendall" = "kendall"
          )
        ),
        selectInput(
          inputId = "mdsDistMethod",
          label = "Distance Measure",
          choices = c(
            "Euclidean" = "euclidean",
            "Maximum" = "maximum",
            "Manhattan" = "manhattan",
            "Canberra" = "canberra",
            "Binary" = "binary",
            "Minkowski" = "minkowski"
          )
        ),
        selectInput(
          inputId = "mds",
          label = "MDS Method",
          choices = c(
            "Classical MDS" = "Classical MDS",
            "Nonmetric MDS" = "Nonmetric MDS"
          )
        )
      ),
      column(9, plotlyOutput("mdsPlotObject") %>% withSpinner())
    ))
  } else {
    helpText("No data for ploting. Please import dataset and assign group information first.")
  }
})

# PCA Plot Scree ----
output$pcaPlotObjectScree <- renderPlotly({
  if (length(variables$count.data) > 0) {
    dt <- variables$count.data
    if (input$pcaTransform == TRUE) {
      data <- log1p(dt)
    } else {
      data <- dt
    }
    data <- data[apply(data, 1, var) != 0, ]
    if (!is.na(input$pcaTopGene) & input$pcaTopGene < nrow(data)) {
      data <- t(data[order(apply(data, 1, var), decreasing = TRUE)[1:input$pcaTopGene], ])
    }
    
    data.pca.all <- prcomp(data,
                           center = input$pcaCenter,
                           scale. = input$pcaScale
    )
    summaryTable <- summary(data.pca.all)$importance
    if(ncol(summaryTable) > 8){
      summaryTable <- summaryTable[,1:8]
    }
    p <- plot_ly(
      x = colnames(summaryTable),
      y = summaryTable[2, ],
      text = paste0(summaryTable[2, ] * 100, "%"),
      textposition = "auto",
      type = "bar",
      name = "Proportion of Variance"
    ) %>%
      add_trace(
        y = summaryTable[3, ],
        type = "scatter",
        mode = "lines+markers",
        name = "Cumulative Proportion"
      ) %>%
      layout(
        xaxis = list(title = "Principal Components"),
        yaxis = list(
          title = "Proportion of Variance",
          tickformat = "%"
        ),
        title = "Scree Plot",
        legend = list(
          orientation = "h",
          xanchor = "center",
          x = 0.5,
          y = 1.05
        )
      ) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "svg",
          filename = "Scree_Plot"
        )
      )
    variables$screePlot <- p
    p
  } else {
    return(0)
  }
})
# PCA Plot 3D ----
output$pcaPlotObject3d <- renderPlotly({
  if (length(variables$count.data) > 0) {
    dt <- variables$count.data
    if (input$pcaTransform == TRUE) {
      data <- log1p(dt)
    } else {
      data <- dt
    }
    data <- data[apply(data, 1, var) != 0, ]
    if (!is.na(input$pcaTopGene) & input$pcaTopGene < nrow(data)) {
      data <- t(data[order(apply(data, 1, var), decreasing = TRUE)[1:input$pcaTopGene], ])
    }
    data.pca.all <- prcomp(data,
                           center = input$pcaCenter,
                           scale. = input$pcaScale
    )
    
    data <- data.frame(data.pca.all$x)
    data$name <- rownames(data)
    group <- variables$group_import
    group$name <- rownames(group)
    data <- left_join(x = data, y = group, by = "name")
    importance <- summary(data.pca.all)$importance
    
    #label sample names or not
    if(input$pcaLabel){
      mode <- "markers+text"
    }else{
      mode <- "markers"
    }
    
    
    p <- plot_ly(
      data = data,
      x = ~PC1,
      y = ~PC2,
      z = ~PC3,
      color = ~ factor(group),
      text = ~name,
      textposition = "top right",
      type = "scatter3d",
      mode = mode
    ) %>%
      layout(title = "PCA Plot (3D)",
             scene = list(
               xaxis = list(title = (paste0("PC1(",round(importance[2,1]*100,digits = 2),"%)"))),
               yaxis = list(title = (paste0("PC2(",round(importance[2,2]*100,digits = 2),"%)"))),
               zaxis = list(title = (paste0("PC3(",round(importance[3,2]*100,digits = 2),"%)")))
             )
             ) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "svg",
          filename = "PCA_Plot_in_3D"
        )
      )
    variables$pca3d <- p
    p
  } else {
    return(0)
  }
})
# PCA Plot 2D ----
output$pcaPlotObject2d <- renderPlotly({
  if (length(variables$count.data) > 0) {
    dt <- variables$count.data
    if (input$pcaTransform == TRUE) {
      data <- log1p(dt)
    } else {
      data <- dt
    }
    data <- data[apply(data, 1, var) != 0, ]
    if (!is.na(input$pcaTopGene) & input$pcaTopGene < nrow(data)) {
      data <- t(data[order(apply(data, 1, var), decreasing = TRUE)[1:input$pcaTopGene], ])
    }
    
   #perform PCA
    data.pca.all <- prcomp(data,
                           center = input$pcaCenter,
                           scale. = input$pcaScale
    )
    data <- data.frame(data.pca.all$x)
    data$name <- rownames(data)
    group <- variables$group_import
    group$name <- rownames(group)
    data <- left_join(x = data, y = group, by = "name")
    importance <- summary(data.pca.all)$importance
    
    
    #label sample names or not
    if(input$pcaLabel){
      mode <- "markers+text"
    }else{
      mode <- "markers"
    }
    
    #plot PCA
    p <- plot_ly(
      data = data,
      x = ~PC1,
      y = ~PC2,
      color = ~ factor(group),
      text = ~name,
      type = "scatter",
      mode = mode,
      marker = list(size = 8)
    ) %>%
      layout(title = "PCA Plot (2D)",
             xaxis = list(title = (paste0("PC1(",round(importance[2,1]*100,digits = 2),"%)"))),
             yaxis = list(title = (paste0("PC2(",round(importance[2,2]*100,digits = 2),"%)")))
      ) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "svg",
          filename = "PCA_Plot2D"
        )
      )

    variables$pca2d <- p
    p
  } else {
    return(0)
  }
})

# PCA Summary Table ----
output$pcaSummaryObject <- DT::renderDataTable({
  if (length(variables$count.data) > 0) {
    dt <- variables$count.data
    if (input$pcaTransform == TRUE) {
      data <- log1p(dt)
    } else {
      data <- dt
    }
    data <- data[apply(data, 1, var) != 0, ]
    if (!is.na(input$pcaTopGene) & input$pcaTopGene < nrow(data)) {
      data <- t(data[order(apply(data, 1, var), decreasing = TRUE)[1:input$pcaTopGene], ])
    }
    data.pca.all <- prcomp(data,
                           center = input$pcaCenter,
                           scale. = input$pcaScale
    )
    
    variables$pcaParameter <- list(
      "pcaTransform" = input$pcaTransform,
      "pcaCenter" = input$pcaCenter,
      "pcaScale" = input$pcaScale,
      "pcaTopGene" = input$pcaTopGene
    )
    
    summaryTable <- summary(data.pca.all)$importance
    row.names(summaryTable)[1] <- "Standard Deviation"
    summaryTable <- t(summaryTable)
    t <- DT::datatable(summaryTable, options = list(
      dom = "Bt",
      buttons = list(
        "copy",
        "print",
        list(
          extend = "collection",
          buttons = c("csv", "excel", "pdf"),
          text = "Download"
        )
      )
    )) %>%
      formatRound(
        columns = colnames(summaryTable),
        digits = 3
      ) %>%
      formatStyle(
        "Proportion of Variance",
        background = styleColorBar(range(0, 1), "lightblue"),
        backgroundSize = "98% 88%",
        backgroundRepeat = "no-repeat",
        backgroundPosition = "center"
      ) %>%
      formatStyle(
        "Standard Deviation",
        background = styleColorBar(range(0, summaryTable[, 1]), "lightblue"),
        backgroundSize = "98% 88%",
        backgroundRepeat = "no-repeat",
        backgroundPosition = "center"
      ) %>%
      formatStyle(
        "Cumulative Proportion",
        background = styleColorBar(range(0, 1), "lightblue"),
        backgroundSize = "98% 88%",
        backgroundRepeat = "no-repeat",
        backgroundPosition = "center"
      )
    variables$summaryPCA <- t
    t
  } else {
    return()
  }
})

output$pcatopProteinPercent <- renderUI({
  dt <- variables$count.data
  selected <- input$pcaTopGene
  dim <- nrow(dt)
  percent <- selected/dim*100
  percent <- round(percent,digits = 2)
  
  tagList(
    tipify(
      tags$p(tags$b("Percentage"), ":", percent,"%"),
      title = "Percentage of total proteins",
      placement = "left"
    ))
  
})

# render PCA UI ----
output$pcaUI <- renderUI({
  if (dataImportCheck$importRunValue) {
    tagList(fluidRow(
      column(
        3,
        helpText(
          "PCA is performed on the top n proteins selected by none-zero row variance (or the most variable proteins)."
        ),
        tipify(
          numericInput(
            inputId = "pcaTopGene",
            label = "Top Variable Proteins",
            value = 100,
            min = -1,
            step = 1
          ),
          title = "How many of the most variable proteins should be used for calculating the PCA. Use all none-zero row variance protein if none value is supplied.",
          placement = "left"
        ),
        uiOutput("pcatopProteinPercent"),
        tipify(
          materialSwitch(
            inputId = "pcaTransform",
            label = "Log2 transform",
            value = TRUE,
            right = TRUE,
            status = "primary"
          ),
          title = "Whether the raw intensity should be performed log2 transformation before analysis.",
          placement = "left"
        ),
        tipify(
          materialSwitch(
            inputId = "pcaCenter",
            label = "Center",
            value = TRUE,
            right = TRUE,
            status = "primary"
          ),
          title = "Whether the value should be shifted to be zero centered.",
          placement = "left"
        ),
        tipify(
          materialSwitch(
            inputId = "pcaScale",
            label = "Scale",
            value = TRUE,
            right = TRUE,
            status = "primary"
          ),
          title = "Whether the value should be scaled to have unit variance before the analysis.",
          placement = "left"
        ),
        tipify(
          materialSwitch(
            inputId = "pcaLabel",
            label = "Show labels",
            value = FALSE,
            right = TRUE,
            status = "primary"
          ),
          title = "Whether sample names should be displayed on plot",
          placement = "left"
        )
      ),
      column(
        9,
        tabsetPanel(
          tabPanel(title = "2D Plot", plotlyOutput("pcaPlotObject2d") %>% withSpinner()),
          tabPanel(title = "3D Plot", plotlyOutput("pcaPlotObject3d") %>% withSpinner()),
          tabPanel(
            title = "Summary Table",
            DT::dataTableOutput("pcaSummaryObject") %>% withSpinner()
          ),
          tabPanel(
            title = "Scree Plot",
            plotlyOutput("pcaPlotObjectScree") %>% withSpinner()
          )

        )
      )
    ))
  } else {
    helpText("No data for ploting. Please import dataset and assign group information first.")
  }
})

# Sample-Sample correlation ----
output$dendPlotObject <- renderPlotly({
  if (length(variables$count.data) > 0) {
    dt <- variables$count.data
    
    dt <- dt[apply(dt, 1, var) != 0, ]
    if (!is.na(input$dendTopGene) & input$dendTopGene < nrow(dt)) {
      data <- dt[order(apply(dt, 1, var), decreasing = TRUE)[1:input$dendTopGene], ]
    }
    
    data <- data[rowSums(data) > 0, ]
    data <- data.frame(cor(data, method = input$dendCor))
    data.cl.count <- length(unique(variables$group_import$group))
    p <- heatmaply(
      data,
      k_col = data.cl.count,
      k_row = data.cl.count,
      hclust_method = input$dendCluster,
      labRow = rownames(data),
      labCol = colnames(data),
      colors = GnBu(100)
    ) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "svg",
          filename = "Hierarchical_Clustering"
        )
      )
    variables$original_heatmap <- p
    p
  } else {
    return()
  }
})

# Render dend and heatmap UI ----

output$dendTopProteinPreview <- renderUI({
  dt <- variables$count.data
  selected <- input$dendTopGene
  dim <- nrow(dt)
  percent <- selected/dim*100
  percent <- round(percent,digits = 2)
  
  tagList(
    tipify(
      tags$p(tags$b("Percentage"), ":", percent,"%"),
      title = "Percentage of total proteins",
      placement = "left"
    ))
  
})
output$dendUI <- renderUI({
  if (dataImportCheck$importRunValue) {
    tagList(fluidRow(
      column(
        3,
        tipify(
          numericInput(
            inputId = "dendTopGene",
            label = "Top Variable Proteins",
            value = 100,
            min = -1,
            step = 1
          ),
          title = "How many of the most variable proteins should be used for calculating correlation. Use all none-zero row variance protein if none value is supplied.",
          placement = "left"
        ),
        uiOutput("dendTopProteinPreview"),
        selectInput(
          inputId = "dendCluster",
          label = "Agglomeration Method",
          choices = c(
            "Complete" = "complete",
            "Ward.D" = "ward.D",
            "Ward.D2" = "ward.D2",
            "Single" = "single",
            "UPGMA (Average)" = "average",
            "WPGMA (Mcquitty)" = "mcquitty" 
          )
        ),
        selectInput(
          inputId = "dendCor",
          label = "Distance Measure",
          choices = c(
            "Spearman" = "spearman",
            "Pearson" = "pearson"
          )
        )
      ),
      column(9, plotlyOutput("dendPlotObject") %>% withSpinner())
    ))
  } else {
    helpText("No data for ploting. Please import dataset and assign group information first.")
  }
})
