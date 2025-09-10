# server-normalization.R

# inputs:
#method
#raw_exp
#batch_vector
#output
#internal_standard # required if (method eq 'InternalStandard')
#----------------------------------------------------------------------------
library(limma)
library(dplyr)
#source('commandLineArgs.R')

# normalize Batche effects By Internal Standard
normalizeBatchesByInternalStandard = function(tb, bch, std) {
  norm = tb
  goldStd = as.character(std[std[, 3] == 'standard', 2])
  for (i in 1:nrow(std)) {
    if (std[i, 3] == 'non_standard') {
      bchStandard = as.character(std[i, 2])
      f = tb[, goldStd] - tb[, bchStandard]
      # print(f)
      f[is.na(f)] <- 0
      # print(f)
      crt = as.character(std[i, 1])
      #smp=grep(crt,bch) # bug here
      smp = which(bch == crt)# corrected on 6/19/19
      norm[, smp] = tb[, smp] + f
    }
  }
  return(norm)
}
#-----------------------------------------------------------------------------
# To parse R parameters
# Credit to Nurcan Tuncbag from MIT
#args=(commandArgs(TRUE))
for (e in commandArgs(T)) {
  ta = strsplit(e, "=", fixed = TRUE)
  var = ta[[1]][1]
  if (!is.na(ta[[1]][2])) {
    temp = ta[[1]][2]
    var = substr(ta[[1]][1], 2, nchar(ta[[1]][1]))
    if (substr(ta[[1]][1], nchar(ta[[1]][1]), nchar(ta[[1]][1])) == "I") {
      temp = as.integer(temp)
    }
    if (substr(ta[[1]][1], nchar(ta[[1]][1]), nchar(ta[[1]][1])) == "N") {
      temp = as.numeric(temp)
    }
    if (substr(ta[[1]][1], nchar(ta[[1]][1]), nchar(ta[[1]][1])) == "V") {
      temp = strsplit(temp, ',')[[1]]
    }
    assign(var, temp)
    cat("assigned ", var, " the value of |", temp, "|\n")
  } else {
    var_fields = strsplit(var, '-')[[1]]
    var = var_fields[length(var_fields)]
    assign(var, TRUE)
    cat("assigned ", var, " the value of TRUE\n")
  }
}

#-----------------------------------------------------------------------------


#perform normalization process

normRun <- reactiveValues(normRunValue = FALSE)
#v$importActionValue <- F

output$batchSelect <- renderUI({
  metaData1 <- variables$group
  df = metaData1
  if (is.null(colnames(df))) {
    vars = NULL
  } else {
    vars = colnames(df)[2:ncol(df)]
  }
  
  selectInput("group_batch",
              "Batch selection",
              choices = vars,
              selected = vars[1])
})

output$internalReference <- renderUI({
  metaData1<- variables$group
  df = metaData1
  if (is.null(colnames(df))) {
    vars = NULL
  } else {
    vars = colnames(df)[2:ncol(df)]
  }
  if (input$normalization == "internalNM" || input$normalization == "internal_linearNM") {
    selectInput(
      inputId = "intRefCol",
      label = "Internal reference column",
      choices = vars,
      selected = vars[1]
    )
  }
})


#if the batchlist is uploaded
observeEvent(input$confirmedBatchList, {
  
  rawdata <- variables$CountData
  
  
  
  #run batch effect elimination
  tb = rawdata
  tb <- tb[c(-1, -2)]
  
  smallValue = 1
  tb = log2(tb + smallValue)
  
  
  
  
  
  #input the batch group list from user
  batchgroup <- input$group_batch
  
  groupList <- variables$group
  
  idx <- as.numeric(grep(batchgroup, colnames(groupList)))
  if (colnames(groupList)[1] == "Sample") {
    colnames(groupList)[1] <- tolower(colnames(groupList)[1])
  }
  
  #get the new group
  variables$batchGroupList <-split(groupList$sample, groupList[, idx])
  #names(variables$batchGroupList) <- unique(groupList[, idx])
  
  #get the batch group
  batchgroup <- groupList[, c(1, idx)]
  colnames(batchgroup) <- c("sample", "batch")
  

  
  
  batch_vector = as.data.frame(batchgroup$batch)
  bch = as.vector(batch_vector[, 1])

  
  
  # Normalize data based on selection
  
  #get the internal information
  if (input$normalization == 'internalNM' || input$normalization=='internal_linearNM') {
    idx <- as.numeric(grep(input$intRefCol, colnames(groupList)))
    interRef <- groupList[, idx]
    trans = as.data.frame(cbind(batchgroup$batch, batchgroup$sample, interRef))
    standard_vector = trans[trans$interRef == "internal", ]
    standard_vector$interRef[1] = "standard"
    standard_vector$interRef[2:nrow(standard_vector)] <- "non_standard"
    
    norm = normalizeBatchesByInternalStandard(tb, bch, standard_vector)
    if(input$normalization == 'internal_linearNM'){
      norm = removeBatchEffect(norm, bch)
    }

    
  } else{
    if (input$normalization == 'linearNM') {
      norm = removeBatchEffect(tb, bch)
    }
  }
  
  norm <- 2 ^ norm
  
  #save the normalized data
  variables$normedData <- cbind(variables$CountData[, c(1, 2)], norm)
  
  data.cl <- rep(0, ncol(variables$normedData))
  for (i in seq_along(variables$batchGroupList)) {
    indices <- match(variables$batchGroupList[[i]], colnames(variables$normedData))
    data.cl[indices] <- names(variables$batchGroupList)[i]
  }
  # for (i in 1:length(variables$batchGroupList)) {
  #   data.cl[unlist(lapply(variables$batchGroupList[[i]], convert2cl, df = variables$normedData))] <- names(variables$batchGroupList[i])
  # }
  
  # Storage convert group list to local
  variables$batchGroupListConvert <- data.cl
  
  
  
  variables$normed.count.data <- variables$normedData[, data.cl != 0]
  variables$normed.count.data <- as.matrix(variables$normed.count.data)
  
  group2 <- data.cl[data.cl != 0]
  group2 <- data.frame(group = group2)
  group2 <- data.frame(lapply(group2, as.factor))
  rownames(group2) <- colnames(variables$count.data)
  variables$normed_group_import <- group2
  
  
  normRun$normRunValue <- TRUE
  
  
  tbtoshow <- round(norm, digits = 2)
  tbtoshow <- cbind(variables$CountData[, c(1, 2)], tbtoshow)
  # Render normalized result table on the right top ----
  output$downloadnorm <- downloadHandler(
    filename = "normalized_data.csv",
    content = function(file) {
      write.csv(variables$normedData, file)
    }
  )
  
  output$normTable <- DT::renderDataTable({
    DT::datatable(
      tbtoshow,
      colnames = c("Accession Number" = 1),
      extensions = c("Scroller", "RowReorder"),
      option = list(
        rowReorder = TRUE,
        deferRender = TRUE,
        scrollY = 400,
        scroller = TRUE,
        scrollX = TRUE,
        searchHighlight = TRUE,
        orderClasses = TRUE
      )
    )
  }, server = T)
  
  
  
  
})

# Render DataTable of row data proteins ----
output$normalizationResultTable <- renderUI({
  if(normRun$normRunValue){
    tagList(
      fluidRow(column(
        12, 
        downloadButton("downloadnorm", "Download Result Table"),
        DT::dataTableOutput('normTable') %>% withSpinner()
      )))} else {
        helpText("Click [Run Normalization] to obtain Result Table.")
      }
})






# Output different plots on bottom right

# This function render a boxplot of sample distribution ----
output$normed_sampleDistributionBox <- renderPlotly({
  if (length(variables$normed.count.data) > 0) {
    normed_data <- variables$normed.count.data
    
    normed_data <- as.data.frame(log2(normed_data))
    
    # Filter
    log2data <- normed_data %>% filter_all(all_vars(. > input$normed_sampleDistributionFilterLow |
                                                      is.na(.)))
    
    data_stack <- data.frame(stack(log2data))
    
    # Add Group info
    group <-
      data.frame(
        "ind" = rownames(variables$normed_group_import),
        "group" = variables$normed_group_import$group
      )
    
    data <- left_join(data_stack, group, by = "ind")
    data <- arrange(data, group)
    
    p <- plot_ly(
      data = data,
      x = ~ ind,
      y = ~ values,
      type = "box",
      split = ~ group,
      color = ~ group
    ) %>%
      layout(
        title = input$norm_sampleDistributionTitle,
        xaxis = list(
          title = input$norm_sampleDistributionXlab,
          categoryarray = "array",
          categoryarray = ~ col
        ),
        yaxis = list(title = input$norm_sampleDistributionYlab)
      ) %>%
      config(
        toImageButtonOptions = list(
          format = "svg",
          filename = input$norm_sampleDistributionTitle
        )
      )
    variables$normed_sampleDistributionBar <- p
    p
  } else {
    return()
  }
})

# Render a density plot of sample distribution ----
output$normed_sampleDistributionDensity <- renderPlotly({
  if (length(variables$normed.count.data) > 0) {
    normed_data <- variables$normed.count.data
    if (input$norm_densityFilter != "Do not filter") {
      count <-
        normed_data[rowSums(normed_data) > as.numeric(input$norm_densityFilter), ]
    } else {
      count <- normed_data
    }
    data <- log2(count + 1)
    
    group <- variables$normed_group_import
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
        title = input$norm_sampleDistributionDenstityTitle,
        xaxis = list(title = input$norm_sampleDistributionDensityXlab),
        yaxis = list(title = input$norm_sampleDistributionDensityYlab)
      ) %>%
      config(
        toImageButtonOptions = list(
          format = "svg",
          filename = input$norm_sampleDistributionDenstityTitle
        )
      )
    variables$normed_sampleDistributionDensity <- pp
    pp
  } else {
    return()
  }
})

# Sample Distribution Density Plot UI ----
output$norm_sampleDistributionDensityPanel <- renderUI({
  if (normRun$normRunValue) {
    tagList(fluidRow(
      column(
        3,
        popify(
          helpText(
            "Filter proteins with a total read count smaller than thresholds."
          ),
          title = "Reference",
          content = 'Sultan, Marc, et al. <a href="http://science.sciencemag.org/content/321/5891/956">"A global view of protein activity and alternative splicing by deep sequencing of the human transcriptome."</a> <i>Science</i> 321.5891 (2008): 956-960.',
          placement = "left"
        ),
        sliderTextInput(
          inputId = "norm_densityFilter",
          label = "Filter proteins threshold",
          choices = c("Do not filter", c(0:30))
        ),
        textInput(
          inputId = "norm_sampleDistributionDenstityTitle",
          label = "Title",
          value = "Raw Intensity",
          placeholder = "Raw Intensity"
        ),
        textInput(
          inputId = "norm_sampleDistributionDensityXlab",
          label = "X label",
          value = "log<sub>2</sub>(Intensity + 1)",
          placeholder = "log<sub>2</sub>(Intensity + 1)"
        ),
        textInput(
          inputId = "norm_sampleDistributionDensityYlab",
          label = "Y label",
          value = "Density",
          placeholder = "Density"
        )
      ),
      column(
        9,
        plotlyOutput("normed_sampleDistributionDensity") %>% withSpinner()
      )
    ))
  } else {
    helpText("No data for ploting. Please import dataset and assign group information first.")
  }
})


# Sample Distribution Boxplot UI ----
output$norm_DistributionBoxPanel <- renderUI({
  if (normRun$normRunValue) {
    tagList(fluidRow(
      column(
        3,
        sliderInput(
          inputId = "normed_sampleDistributionFilterLow",
          label = "Filter low proteins",
          min = 0,
          max = 30,
          value = 0
        ),
        textInput(
          inputId = "sampleDistributionTitle",
          label = "Title",
          value = "Original Raw Intensity",
          placeholder = "Original Raw Intensity"
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
          value = "log<sub>2</sub>(Intensity + 1)",
          placeholder = "log<sub>2</sub>(Intensity + 1)"
        )
      ),
      column(
        9,
        plotlyOutput("normed_sampleDistributionBox") %>% withSpinner()
      )
    ))
  } else {
    helpText("No data for ploting. Please import dataset and assign group information first.")
  }
})


# Filtering low count under different low number, barplot ----
output$norm_lowCountFilterByCutoff <- renderPlotly({
  if (length(variables$normed.count.data) > 0) {
    normed_data <- variables$normed.count.data
    originalCount <- nrow(normed_data)
    
    lowCount <-
      sapply(0:input$lowCountSlide, function(x) {
        sum(rowSums(normed_data) > x)
      })
    
    lowCountdt <- data.frame(
      "Cutoff" = 0:input$lowCountSlide,
      "Filtered" = originalCount - lowCount,
      "Remain" = lowCount
    )
    
    plot_ly(
      lowCountdt,
      name = "Remain",
      x = ~ Cutoff,
      y = ~ Remain,
      text = ~ Remain,
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
        y = ~ Filtered,
        text = ~ Filtered,
        type = "bar",
        textposition = "inside"
      ) %>%
      layout(
        title = "Filtering Threshold for Low Intensity",
        barmode = "stack",
        xaxis = list(title = "Filtering Low Intensity Cut off"),
        yaxis = list(title = "Protein number")
      ) %>%
      config(
        toImageButtonOptions = list(format = "svg", filename = "Filtering_Threshold_for_Low_Count_Proteins")
      )
  } else {
    return()
  }
})

# Render Filtering Cutoff UI ----
output$norm_lowCountFilterByCutoffUI <- renderUI({
  if (normRun$normRunValue) {
    tagList(fluidRow(
      column(
        3,
        popify(
          helpText(
            "Filter proteins with tota intensity smaller than thresholds."
          ),
          placement = "left"
        ),
        sliderInput(
          inputId = "lowCountSlide",
          label = "Max threshold",
          min = 1,
          max = 35,
          value = 20,
          step = 1
        )
      ),
      column(
        9,
        plotlyOutput("norm_lowCountFilterByCutoff") %>% withSpinner()
      )
    ))
  } else {
    helpText("No data for ploting. Please import dataset and assign group information first.")
  }
})

# PCA Plot Scree ----
output$normed_pcaPlotObjectScree <- renderPlotly({
  if (length(variables$normed.count.data) > 0) {
    normed_data <- variables$normed.count.data
    if (input$normed_pcaTransform == TRUE) {
      normed_data <- log1p(normed_data)
    } else {
      normed_data <- normed_data
    }
    normed_data <- normed_data[apply(normed_data, 1, var) != 0, ]
    if (!is.na(input$normed_pcaTopprotein) &
        input$normed_pcaTopprotein < nrow(normed_data)) {
      normed_data <- t(normed_data[order(apply(normed_data, 1, var), decreasing = TRUE)[1:input$normed_pcaTopprotein], ])
    }
    
    data.pca.all <- prcomp(
      normed_data,
      center = input$normed_pcaCenter,
      scale. = input$normed_pcaScale
    )
    summaryTable <- summary(data.pca.all)$importance
    if (ncol(summaryTable) > 8) {
      summaryTable <- summaryTable[, 1:8]
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
        yaxis = list(title = "Proportion of Variance", tickformat = "%"),
        title = "Scree Plot",
        legend = list(
          orientation = "h",
          xanchor = "center",
          x = 0.5,
          y = 1.05
        )
      ) %>%
      config(toImageButtonOptions = list(format = "svg", filename = "Scree_Plot"))
    variables$normed_screePlot <- p
    p
  } else {
    return(0)
  }
})

# PCA Plot 3D ----
output$normed_pcaPlotObject3d <- renderPlotly({
  if (length(variables$normed.count.data) > 0) {
    data <- variables$normed.count.data
    if (input$normed_pcaTransform == TRUE) {
      data <- log1p(data)
    } else {
      data <- data
    }
    data <- data[apply(data, 1, var) != 0, ]
    if (!is.na(input$normed_pcaTopprotein) &
        input$normed_pcaTopprotein < nrow(data)) {
      data <- t(data[order(apply(data, 1, var), decreasing = TRUE)[1:input$normed_pcaTopprotein], ])
    }
    data.pca.all <- prcomp(data,
                           center = input$normed_pcaCenter,
                           scale. = input$normed_pcaScale)
    
    data <- data.frame(data.pca.all$x)
    data$name <- rownames(data)
    group <- variables$normed_group_import
    group$name <- rownames(group)
    data <- left_join(x = data, y = group, by = "name")
    importance <- summary(data.pca.all)$importance
    
    #label sample names or not
    if (input$normed_pcaLabel) {
      mode <- "markers+text"
    } else{
      mode <- "markers"
    }
    
    p <- plot_ly(
      data = data,
      x = ~ PC1,
      y = ~ PC2,
      z = ~ PC3,
      color = ~ factor(group),
      text = ~ name,
      textposition = "top right",
      type = "scatter3d",
      mode = mode
    ) %>%
      layout(title = "PCA Plot (3D)",
             scene = list(
               xaxis = list(title = (paste0(
                 "PC1(", round(importance[2, 1] * 100, digits = 2), "%)"
               ))),
               yaxis = list(title = (paste0(
                 "PC2(", round(importance[2, 2] * 100, digits = 2), "%)"
               ))),
               zaxis = list(title = (paste0(
                 "PC3(", round(importance[3, 2] * 100, digits = 2), "%)"
               )))
             )) %>%
      config(toImageButtonOptions = list(format = "svg", filename = "PCA_Plot_in_3D"))
    variables$normed_pca3d <- p
    p
  } else {
    return(0)
  }
})
# PCA Plot 2D ----
output$normed_pcaPlotObject2d <- renderPlotly({
  if (length(variables$normed.count.data) > 0) {
    data <- variables$normed.count.data
    if (input$normed_pcaTransform == TRUE) {
      data <- log1p(data)
    } else {
      data <- data
    }
    data <- data[apply(data, 1, var) != 0, ]
    if (!is.na(input$normed_pcaTopprotein) &
        input$normed_pcaTopprotein < nrow(data)) {
      data <- t(data[order(apply(data, 1, var), decreasing = TRUE)[1:input$normed_pcaTopprotein], ])
    }
    data.pca.all <- prcomp(data,
                           center = input$normed_pcaCenter,
                           scale. = input$normed_pcaScale)
    data <- data.frame(data.pca.all$x)
    data$name <- rownames(data)
    group <- variables$normed_group_import
    group$name <- rownames(group)
    data <- left_join(x = data, y = group, by = "name")
    importance <- summary(data.pca.all)$importance
    
    
    #label sample names or not
    if (input$normed_pcaLabel) {
      mode <- "markers+text"
    } else{
      mode <- "markers"
    }
    
    p <- plot_ly(
      data = data,
      x = ~ PC1,
      y = ~ PC2,
      color = ~ factor(group),
      text = ~ name,
      textposition = "top right",
      type = "scatter",
      mode = mode,
      marker = list(size = 8)
    ) %>%
      layout(
        title = "PCA Plot (2D)",
        xaxis = list(title = (paste0(
          "PC1(", round(importance[2, 1] * 100, digits = 2), "%)"
        ))),
        yaxis = list(title = (paste0(
          "PC2(", round(importance[2, 2] * 100, digits = 2), "%)"
        )))
      ) %>%
      config(toImageButtonOptions = list(format = "svg", filename = "PCA_Plot_in_2D"))
    variables$normed_pca2d <- p
    p
  } else {
    return(0)
  }
})

# PCA Summary Table ----
output$normed_pcaSummaryObject <- DT::renderDataTable({
  if (length(variables$normed.count.data) > 0) {
    data <- variables$normed.count.data
    if (input$normed_pcaTransform == TRUE) {
      data <- log1p(data)
    } else {
      data <- data
    }
    data <- data[apply(data, 1, var) != 0, ]
    if (!is.na(input$normed_pcaTopprotein) &
        input$normed_pcaTopprotein < nrow(data)) {
      data <- t(data[order(apply(data, 1, var), decreasing = TRUE)[1:input$normed_pcaTopprotein], ])
    }
    data.pca.all <- prcomp(data,
                           center = input$normed_pcaCenter,
                           scale. = input$normed_pcaScale)
    
    summaryTable <- summary(data.pca.all)$importance
    row.names(summaryTable)[1] <- "Standard Deviation"
    summaryTable <- t(summaryTable)
    t <- DT::datatable(summaryTable, options = list(
      dom = "Bt",
      buttons = list(
        list(
          extend = "collection",
          buttons = c("csv", "excel"),
          text = "Download"
        )
      )
    )) %>%
      formatRound(columns = colnames(summaryTable), digits = 3) %>%
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

# render PCA UI ----
output$normed_pcaTopproteinPreview <- renderUI({
  dt <- variables$normed.count.data
  selected <- input$normed_pcaTopprotein
  dim <- nrow(dt)
  percent <- selected / dim * 100
  percent <- round(percent, digits = 2)
  
  tagList(tipify(
    tags$p(tags$b("Percentage"), ":", percent, "%"),
    title = "Percentage of total proteins",
    placement = "left"
  ))
  
})
output$norm_pcaUI <- renderUI({
  if (normRun$normRunValue) {
    tagList(fluidRow(
      column(
        3,
        helpText(
          "PCA is performed on the top n varaible proteins without missing values."
        ),
        tipify(
          numericInput(
            inputId = "normed_pcaTopprotein",
            label = "Top protein",
            value = 100,
            min = -1,
            step = 1
          ),
          title = "Number of top varied proteins for calculation",
          placement = "left"
        ),
        uiOutput("normed_pcaTopproteinPreview"),
        tipify(
          materialSwitch(
            inputId = "normed_pcaTransform",
            label = "Log2 transform",
            value = TRUE,
            right = TRUE,
            status = "primary"
          ),
          title = "Use log2 transformed data or original expression data.",
          placement = "left"
        ),
        tipify(
          materialSwitch(
            inputId = "normed_pcaCenter",
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
            inputId = "normed_pcaScale",
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
            inputId = "normed_pcaLabel",
            label = "Show labels",
            value = FALSE,
            right = TRUE,
            status = "primary"
          ),
          title = "Whether sample names should be displayed on plot",
          placement = "left"
        )
      ),
      column(9, tabsetPanel(
        tabPanel(
          title = "2D Plot",
          plotlyOutput("normed_pcaPlotObject2d") %>% withSpinner()
        ),
        tabPanel(
          title = "3D Plot",
          plotlyOutput("normed_pcaPlotObject3d") %>% withSpinner()
        ),
        tabPanel(
          title = "Summary Table",
          DT::dataTableOutput("normed_pcaSummaryObject") %>% withSpinner()
        ),
        tabPanel(
          title = "Scree Plot",
          plotlyOutput("normed_pcaPlotObjectScree") %>% withSpinner()
        )
      ))
    ))
  } else {
    helpText("No data for ploting. Please import dataset and assign group information first.")
  }
})

# sample-sample correlation heatmap ----
output$normed_dendPlotObject <- renderPlotly({
  if (length(variables$normed.count.data) > 0) {
    data <- variables$normed.count.data
    
    data <- data[apply(data, 1, var) != 0, ]
    if (!is.na(input$normed_dendTopGene) &
        input$normed_dendTopGene < nrow(data)) {
      data <- data[order(apply(data, 1, var), decreasing = TRUE)[1:input$normed_dendTopGene], ]
    }
    
    data <- data[rowSums(data) > 0, ]
    data <- data.frame(cor(data, method = input$normed_dendCor))
    data.cl.count <- length(unique(variables$normed_group_import$group))
    p <- heatmaply(
      data,
      k_col = data.cl.count,
      k_row = data.cl.count,
      hclust_method = input$normed_dendCluster,
      labRow = rownames(data),
      labCol = colnames(data),
      colors = rev(RdBu(100))
    ) %>%
      config(toImageButtonOptions = list(format = "svg", filename = "Sample_Correlation"))
    variables$normed_heatmap <- p
    p
  } else {
    return()
  }
})

# Render dend and heatmap UI ----
output$normed_dendtopProteinPreview <- renderUI({
  dt <- variables$normed.count.data
  selected <- input$normed_dendTopGene
  dim <- nrow(dt)
  percent <- selected / dim * 100
  percent <- round(percent, digits = 2)
  
  tagList(tipify(
    tags$p(tags$b("Percentage"), ":", percent, "%"),
    title = "Percentage of total proteins",
    placement = "left"
  ))
})

output$norm_dendUI <- renderUI({
  if (normRun$normRunValue) {
    tagList(fluidRow(
      column(
        3,
        tipify(
          numericInput(
            inputId = "normed_dendTopGene",
            label = "Top Variable Proteins",
            value = 100,
            min = -1,
            step = 1
          ),
          title = "How many of the most variable proteins should be used for calculating correlation. Use all none-zero row variance protein if none value is supplied.",
          placement = "left"
        ),
        uiOutput("normed_dendtopProteinPreview"),
        selectInput(
          inputId = "normed_dendCluster",
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
          inputId = "normed_dendCor",
          label = "Distance Measure",
          choices = c("Spearman" = "spearman", "Pearson" = "pearson")
        )
      ),
      column(
        9,
        plotlyOutput("normed_dendPlotObject") %>% withSpinner()
      )
    ))
  } else {
    helpText("No data for ploting. Please import dataset and assign group information first.")
  }
})
