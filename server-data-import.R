# server-data-import.R

options(shiny.maxRequestSize = 60*1024^2)

dataImportCheck <- reactiveValues(importRunValue = FALSE)


# If Upload Data button has been clicked, load the data via upload ----

#validate if there are any duplicated uniprot ids
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

observeEvent(input$uploadExpressionData, {
  showNotification("Uploading file...", type = "message")
  
  shinyCatch(
    {
      #read in uploaded data
      rawdata <-
        data.frame(fread(input$uploadExpressionData$datapath))
      
      #omit NA values
      nNA <- sum(rowSums(is.na(rawdata))>0)
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
        first <- batchidx[length(batchidx)]+5
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
  #validate(dupValidate(new))
  new <- new[!duplicated(new$pid),]
  variables$dup <- nrow(new[duplicated(new$pid),])
  row.names(new) <- new$pid
  new <- new[,-1]
  
  
  #save the raw data as a global variable
  if(input$intensityType  == "log2"){
    data <- new[,3:ncol(new)]
    data_trans <- 2^data
    temp <- cbind(new[,1:2], data_trans)
    new <- temp
  }
  
  variables$CountData <- new
  
  dataImportCheck$importRunValue <- FALSE
  
  showNotification("Data uploaded successfully.", type = "message")
  
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
  
  
  DT::datatable(
    df,
    colnames = c("Accession Number" = 1),
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
  if (nrow(variables$CountData) == 0) {
    tags$p("No data to show. Download", tags$code("Example"), "or", tags$code("Upload"), "your own dataset.")
  } else {
    tagList(
      fluidRow(
        column(
          12,
          DT::dataTableOutput("table"),
          plotlyOutput('raw_boxplot')%>% withSpinner()
        )
      )
    )
    
  }
})

output$raw_boxplot <- renderPlotly({
  req(nrow(variables$CountData) > 0)     # require expression data
  req(input$table_rows_selected)         # require a row click
  
  # Expression data (skip first 2 meta columns)
  data <- variables$CountData[, 3:ncol(variables$CountData)]
  
  # Row index when clicking the row
  rowInd <- input$table_rows_selected
  expr <- as.numeric(data[rowInd, ])
  
  # --- Group mapping logic ---
  if (is.null(variables$group) || is.null(input$groups1)) {
    # No group info uploaded yet → default group
    data.cl <- rep("Samples", ncol(data))
  } else {
    # Use uploaded group info
    data.cl <- variables$groupListConvert
    data.cl <- data.cl[data.cl!=0]
    if (all(data.cl == 0)) {
      data.cl <- rep("Unassigned", ncol(data))
    }
  }
  
  # Build plotting dataframe
  #group_levels <- unique(names(variables$groupList))
  df <- data.frame(
    samples   = colnames(data),
    intensity = round(log2(expr), digits = 2),
    group     = factor(data.cl, levels = unique(data.cl))
  )
  
  # Order x-axis by group
  xOrderVector <- unique(df$samples[order(df$group)])
  xform <- list(
    categoryorder = "array",
    categoryarray = xOrderVector,
    title = ""
  )
  
  # Plot barplot
  plot_ly(
    data = df,
    x = ~samples,
    y = ~intensity,
    color = ~group,          # color by group (or "All Samples")
    text = ~intensity,
    type = "bar",
    showlegend = TRUE
  ) %>%
    layout(
      xaxis = xform,
      yaxis = list(
        title = "Intensity distribution of original data",
        range = c(0.9 * min(df$intensity), max(df$intensity) + 0.5)
      ),
      title = rownames(data)[rowInd]
    ) %>%
    plotly::config(
      toImageButtonOptions = list(
        format = "svg",
        filename = rownames(data)[rowInd]
      )
    )
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
      showNotification("Uploading file...", type = "message")
      
      
      group <- data.frame(fread(input$uploadGroup$datapath, header = T, blank.lines.skip = T))
      group_no_header <- data.frame(fread(input$uploadGroup$datapath,header = FALSE))
      
      showNotification("Data uploaded successfully.", type = "message")
      
      #store group information
      colnames(group)[1] <- "sample"
      variables$group <- group
      
      
      # variables$group_no_header <- group_no_header
      group_no_header <- group_no_header[-1,]
      groupName <- colnames(group)
      
      
      # Create groupList
      variables$groupList <- split(group_no_header$V1, group_no_header$V2)
      
      # Match the groups in group info to expression data samples
      data.cl <- rep(0, ncol(variables$CountData))
      for (i in seq_along(variables$groupList)) {
        indices <- match(variables$groupList[[i]], colnames(variables$CountData))
        data.cl[indices] <- names(variables$groupList)[i]
      }
      
      
      
      # Storage convert group list to local
      variables$groupListConvert <- data.cl
      
      #store numeric data
      variables$count.data <- as.matrix(variables$CountData[, data.cl != 0])
      
      # Check if sample info matches raw data
      if (ncol(variables$count.data) != nrow(group)) {
        sendSweetAlert(
          session = session,
          title = "Warning",
          text = "Sample number should match group info.",
          type = "warning"
        )
        return()
      }
      
      # Store converted group info
      group2 <- data.frame(group = factor(data.cl[data.cl != 0]))
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
  # Import the group selection
  groupList <- as.data.frame(variables$group)
  group <- input$groups1
  
  # Find the column index for the selected group
  idx <- grep(group, colnames(groupList))
  
  # Unify sample column name
  colnames(groupList)[1] <- "sample"
  
  # Create the second group list based on the selected group column
  variables$groupList <- split(groupList$sample, groupList[, idx])
  
  # Match samples in the group with sample in the data
  data.cl <- rep(0, ncol(variables$CountData))
  for (i in seq_along(variables$groupList)) {
    indices <- match(variables$groupList[[i]], colnames(variables$CountData))
    data.cl[indices] <- names(variables$groupList)[i]
  }
  
  # Update converted group list
  variables$groupListConvert <- data.cl
  
  # store only numeric data
  variables$count.data <- as.matrix(variables$CountData[, data.cl != 0])
  
  # Update group info to selected column
  group2 <- data.frame(group = factor(data.cl[data.cl != 0]))
  rownames(group2) <- colnames(variables$count.data)
  variables$group_import <- group2
  
})

#add box plot

output$importDataSummary <- renderUI({
  dt <- datasetInput()
  dup <- variables$dup
  
  rowCount <- nrow(dt)
  groupCount <- length(variables$groupList)
  groupText <- sapply(variables$groupList, length)
  if (length(groupText) > 0) {
    gText <- paste0(names(groupText), ": ", groupText, collapse = "<br>")
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
      tags$p(tags$b("N", tags$sub("Protein/Peptide")), ":", rowCount),
      title = "Number of proteins/peptides",
      placement = "left"
    ),
    tipify(
      tags$p(tags$b("N", tags$sub("group")),": ", groupCount),
      title = "Number of groups",
      placement = "left",
    ),
    tipify(
      tags$p(tags$b("N", tags$sub("samples in group")), HTML(": <br>"), HTML(gText)),
      title = "Number of samples in each group",
      placement = "left"
    ),
    tipify(
      tags$p(tags$b("N", tags$sub("replicates")), ": ", dup),
      title = "Number of replicated gene names",
      placement = "left"
    ),
    tipify(
      tags$p(tags$b("Nrow",tags$sub("NA")), ": ", nNA),
      title = "Rows of proteins/petides that contain at least 1 NA",
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
    
    count <- dt
    data <- log2(count)
    
    group <- variables$group_import
    
    densityTable <- apply(data, 2, function(x) {
      density(na.omit(x))
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
        
        textInput(
          inputId = "sampleDistributionDenstityTitle",
          label = "Title",
          value = "Sample density distribution",
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
          helpText("Filter proteins with intensities smaller than thresholds."),
          title = "Reference",
          content = "Filter proteins by intensity values",
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


# PCA Plot Scree ----
output$pcaPlotObjectScree <- renderPlotly({
  if (length(variables$count.data) > 0) {
    dt <- variables$count.data
    if (input$pcaTransform == TRUE) {
      data <- log2(dt)
    } else {
      data <- dt
    }
    data <- data[apply(data, 1, var) != 0, ]
    if (!is.na(input$pcaTopGene) & input$pcaTopGene < nrow(data)) {
      data <- t(data[order(apply(data, 1, var), decreasing = TRUE)[1:input$pcaTopGene], ])
    }
    data <- data[,colSums(is.na(data))<nrow(data)]
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
      #text = paste0(summaryTable[2, ] * 100, "%"),
      text = sprintf("%.2f%%", summaryTable[2, ] * 100),
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
          tickformat = ".2%"  # Limits to two decimal places
        ),
        title = "Scree Plot",
        legend = list(
          orientation = "h",
          xanchor = "center"
          # x = 0.5,
          # y = 1.05
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
      data <- log2(dt)
    } else {
      data <- dt
    }
    data <- data[apply(data, 1, var) != 0, ]
    if (!is.na(input$pcaTopGene) & input$pcaTopGene < nrow(data)) {
      data <- t(data[order(apply(data, 1, var), decreasing = TRUE)[1:input$pcaTopGene], ])
    }
    data <- data[,colSums(is.na(data))<nrow(data)]
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
      data <- log2(dt)
    } else {
      data <- dt
    }
    data <- data[apply(data, 1, var) != 0, ]
    if (!is.na(input$pcaTopGene) & input$pcaTopGene < nrow(data)) {
      data <- t(data[order(apply(data, 1, var), decreasing = TRUE)[1:input$pcaTopGene], ])
    }
    
    data <- data[,colSums(is.na(data))<nrow(data)]
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
      marker = list(size = 10)
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
      data <- log2(dt)
    } else {
      data <- dt
    }
    data <- data[apply(data, 1, var) != 0, ]
    if (!is.na(input$pcaTopGene) & input$pcaTopGene < nrow(data)) {
      data <- t(data[order(apply(data, 1, var), decreasing = TRUE)[1:input$pcaTopGene], ])
    }
    data <- data[,colSums(is.na(data))<nrow(data)]
    data.pca.all <- prcomp(data,
                           center = input$pcaCenter,
                           scale. = input$pcaScale
    )
    
    
    
    summaryTable <- summary(data.pca.all)$importance
    row.names(summaryTable)[1] <- "Standard Deviation"
    summaryTable <- t(summaryTable)
    displayTable <- DT::datatable(summaryTable,
                                  option = list(
                                    rowReorder = TRUE,
                                    deferRender = TRUE,
                                    autoWidth = TRUE,
                                    pageLength = 20,
                                    searchHighlight = TRUE,
                                    orderClasses = TRUE
                                  )
    ) %>%
      formatRound(
        columns = colnames(summaryTable),
        digits = 3
      ) %>%
      formatStyle(
        "Proportion of Variance",
        background = styleColorBar(range(0, 1), "steelblue"),
        backgroundSize = "100% 90%",
        backgroundRepeat = "no-repeat",
        backgroundPosition = "center"
      ) %>%
      formatStyle(
        "Standard Deviation",
        background = styleColorBar(range(0, summaryTable[, 1]), "steelblue"),
        backgroundSize = "100% 90%",
        backgroundRepeat = "no-repeat",
        backgroundPosition = "center"
      ) %>%
      formatStyle(
        "Cumulative Proportion",
        background = styleColorBar(range(0, 1), "steelblue"),
        backgroundSize = "100% 90%",
        backgroundRepeat = "no-repeat",
        backgroundPosition = "center"
      )
    variables$summaryPCA <- displayTable
    displayTable
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
          "PCA is performed on the top n varaible proteins without missing values."
        ),
        tipify(
          numericInput(
            inputId = "pcaTopGene",
            label = "Top Variable Proteins",
            value = 100,
            min = -1,
            step = 1
          ),
          title = "Number of top varied proteins for calculation",
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
          title = "Use log2 transformed data or original expression data",
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
          title = "Whether values should be zero centered.",
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
          title = "Whether the value should be scaled.",
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
    ngroups <- length(unique(variables$group_import$group))
    p <- heatmaply(
      data,
      k_col = ngroups,
      k_row = ngroups,
      hclust_method = input$dendCluster,
      labRow = rownames(data),
      labCol = colnames(data),
      distfun = input$dendCor,
      colors = rev(RdBu(100))
    ) %>%
      plotly::config(
        toImageButtonOptions = list(
          format = "svg",
          filename = "Sample_Correlation"
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
          title = "Number of most variable proteins used",
          placement = "left"
        ),
        uiOutput("dendTopProteinPreview"),
        selectInput(
          inputId = "dendCluster",
          label = "Clustering method",
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
          label = "Distance function",
          choices = c(
            "Spearman" = "spearman",
            "Pearson" = "pearson"
          )
        )
      ),
      column(9, plotlyOutput("dendPlotObject",height = 700) %>% withSpinner())
    ))
  } else {
    helpText("No data for ploting. Please import dataset and assign group information first.")
  }
})
