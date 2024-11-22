#server-covariance.R

#libraries
library(tidyr)
library(tidyverse)
library(dplyr)
library(stringr)
library(ggplot2)
library(plotly)

#if run correction button is clicked, run the test
corRun <- reactiveValues(corRunValue = FALSE)

output$covVars <- renderUI({
  
  df <- variables$group
  if (is.null(colnames(df))) {
    vars = NULL
  } else {
    vars = colnames(df)[2: ncol(df)]    
  }
  
  selectInput("covGroup", "Covariance group",
              choices = vars,multiple = TRUE)
})

# observeEvent(input$confirmedCovList,{
#   # if (input$CovInput == "") {
#   #   sendSweetAlert(
#   #     session = session,
#   #     title = "ERROR",
#   #     text = "Please input group information!",
#   #     type = "error"
#   #   )
#   #   return()
#   # }
#   # 
#   tryCatch(
#     {
#       progressSweetAlert(
#         session = session,
#         id = "covinput",
#         title = "Retrieving group information...",
#         display_pct = TRUE,
#         value = 0
#       )
#       
# 
#       
#       cov_group_data <- fread(input$CovInput, header = TRUE)
#       #store group information
#       if(colnames(cov_group_data) [1] == "Sample"){
#         colnames(cov_group_data)[1] <- tolower(colnames(cov_group_data)[1])
#       }
#       variables$cov_group <- cov_group_data
#       cov_group <- colnames(cov_group_data)
#       cov_group <- cov_group[,-1]
#     
#       
#     },
#     
#     error = function(e) {
#       sendSweetAlert(
#         session = session,
#         title = "ERROR",
#         text = "Check your group information format!",
#         type = "error"
#       )
#       return()
#     },
#     warning = function(w) {
#       sendSweetAlert(
#         session = session,
#         title = "Group error!",
#         text = "Check your group information format!",
#         type = "error"
#       )
#       return()
#     }
#   )
#   
#   
#   updateProgressBar(
#     session = session,
#     id = "covinput",
#     title = "Loading covriates...",
#     value = 40
#   )
#   
#   # output$covVars <- renderUI({
#   #   df = variables$cov_group
#   #   if (is.null(colnames(df))) {
#   #     Vars = NULL
#   #   } else {
#   #     Vars = colnames(df)[2: ncol(df)]    
#   #   }
#   #   checkboxGroupButtons("selectVars", "Select covariate",
#   #               choices = Vars)
#   # })
#   
#   updateProgressBar(
#     session = session,
#     id = "covinput",
#     title = "Data Processing",
#     value = 60
#   )
#   
#   updateProgressBar(
#     session = session,
#     id = "covinput",
#     title = "Data Processing",
#     value = 100
#   )
#   
#   closeSweetAlert(session = session)
#   sendSweetAlert(
#     session = session,
#     title = "DONE",
#     text = "Covariates uploaded successfully",
#     type = "success"
#   )
#   
# })

observeEvent(input$runCorrection,{
  progressSweetAlert(
    session = session,
    id = "correction",
    title = "Starting correction",
    display_pct = TRUE,
    value = 0
  )

  #select covariates columns
  choice <- input$covGroup
  group <- as.data.frame(variables$group)
  
  if(colnames(group)[1] == "Sample"){
    colnames(group)[1] <- tolower(colnames(group)[1])
  }
  
  #save covariates group information
  group <- group[c("sample",choice)]
  variables$cov_group <- group

  
  #based on inport load normalized data or original data if normalization not required
  if(is.data.frame(variables$normedData) &&nrow(variables$normedData)!=0){
    dt <- variables$normedData
  }else{
    dt <- variables$CountData
  }
  
  dt <- dt[,3:ncol(dt)]
  
  ##check if there are NA values in sample
  row.names(group) <- group$sample
  group_completed <- group[complete.cases(group),]
  group_na <- group[!complete.cases(group),]
  dt_na <- dt[,row.names(group_na)]
  dt <- dt[,row.names(group_completed)]
  
  
  updateProgressBar(
    session = session,
    id = "correction",
    title = "Loading required covriates...",
    value = 30
  )
  
  #make all groups as factors
  index <- names(Filter(is.character,group_completed))
  index <- index[-1]
  group_completed[index] <- lapply(group_completed[index],as.factor)
  #group[index] <- lapply(group[index], function(x) if (is.factor(x)) as.numeric(x) - 1 else x)
  
  updateProgressBar(
    session = session,
    id = "correction",
    title = "Converting non-numeric covriates...",
    value = 40
  )
  
  #get log2 intensities
  intensitylog2 <- log2(dt)
  intensitylog21 <- as.matrix(intensitylog2)
  
  intensitylog2_t <- as.data.frame(t(intensitylog2))
  #
  intensitylog2_t$sample <- rownames(intensitylog2_t)
  #

  
  #
  merged_exp <- merge(group_completed,intensitylog2_t,by = "sample")
  n <- ncol(group_completed)
  n <- n+1
  data_new <- as.data.frame(merged_exp)
  merged_exp <- as.data.frame(merged_exp)
  
  original_intensity <- as.data.frame(merged_exp[,c(1,n:ncol(merged_exp))])
  row.names(original_intensity) <- original_intensity$sample
  original_intensity <- as.data.frame(original_intensity[-1])
  write.csv(merged_exp,"test/merged_exp_test.csv")
  
  updateProgressBar(
    session = session,
    id = "correction",
    title = "Processing data...",
    value = 50
  )
  
  
  #make linear regression expression formula
  frm <- "Expression ~"
  
  if(length(choice)>=2){
    frm <- paste(frm,choice[1])
     # valid.names <- names(choice)[names(choice) != "sample"]  # all but sample
    for(i in 2:length(choice)) {
      
      frm <- paste(frm, "+", choice[i])
    }
    frm <- as.formula(frm)
  }else{
    # valid.names <- names(choice)[names(choice) != "sample"]
    frm <- as.formula(paste("Expression ~",choice[1]))
  }
  
  nstart <- ncol(group_completed)+1
  
  long_data <- pivot_longer(merged_exp,
                            cols = nstart:ncol(merged_exp),
                            names_to = "Protein_ID",
                            values_to = "Expression")

  
  updateProgressBar(
    session = session,
    id = "correction",
    title = "Run Regression...",
    value = 60
  )
  #############################################################
  # Perform linear regression (all at once)
  ######################################################
  
  
  data_new <- merged_exp
  
  # Iterate over each protein ID to fit models and store p-values
  for (protein_id in unique(long_data$Protein_ID)) {
    # Subset data for the current protein
    showNotification(protein_id,type = "message")
    protein_data <- filter(long_data, Protein_ID == protein_id)
    # Fit the model
    model <- shinyCatch(MASS::rlm(frm, data = protein_data,na.action = na.exclude))
    # Check if the model terms exist and extract p-values
    model_summary <- summary(model)
    
    intercept <- model$coefficients[1]
    
    # Get residuals and fit a null model for adjusted data
    residuals_data <- residuals(model)
    protein_data$Expression <- residuals_data+intercept
    
    residuals_data <- as.numeric(residuals_data)
    idx <- grep(protein_id,colnames(merged_exp))
    
    #put corrected data to original data
    data_new[,idx] <- protein_data$Expression
    
  }
  
  updateProgressBar(
    session = session,
    id = "correction",
    title = "Storing data...",
    value = 70
  )
  
  group <- data_new[,1:n-1]
  data_t <- as.data.frame(t(data_new[,c(n:ncol(data_new))]))
  colnames(data_t) <- data_new$sample
  data_f <- 2^data_t
  data_f <- cbind(data_f, dt_na)
  
  variables$corrected_data <- cbind(variables$CountData[,c(1,2)],data_f)
  write.csv(x = data_f,file = "test/after_cor.csv")
  
  # Render result table on the right top ----
  tbtoshow <- round(data_f,digits = 0)
  tbtoshow <- cbind(variables$CountData[,c(1,2)],data_f)
  
  output$cor_resultTable <- DT::renderDataTable({
    if (nrow(tbtoshow) == 0) {
      DT::datatable(tbtoshow)
    } else {
      data = tbtoshow
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
        #   )
          # tags$li(
          #   HTML("<font color=\"#B22222\"><b>Protein ID</b></font> is colored according to FDR cut-off.")
          # )
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
  },server = FALSE)
  
  output$downloadCovResult <- downloadHandler(
    filename = "CovarianceCorrection.csv",
    content = function(file) {
      write.csv(variables$corrected_data, file)
    }
  )
  
  
  updateProgressBar(
    session = session,
    id = "correction",
    title = "Plot Graph",
    value = 80
  )
  
  ##################Plot density Graph###############
  
   #prepare dataset
  nsample <- nrow(group)
  before <- as.data.frame(original_intensity)
  after <- as.data.frame(data_new[,c(n:ncol(data_new))])
  rownames(after) <- data_new$sample
  
  #####################################################
   n_proteins <- dim(intensitylog21)[1]

  
  #prepare dataset to store the correlation and correlation p values
   cor_before <- matrix(nrow = n_proteins, ncol = ncol(group)-1)
   rownames(cor_before) <- rownames(intensitylog2)
   colnames(cor_before) <- colnames(group)[2:ncol(group)]
   
   cor_after <- cor_before
   cor_p_before <- cor_before
   cor_p_after <- cor_before
  
   rownames(group)<-group$sample
   group <- group[-1]
   group[index] <- lapply(group[index], function(x) if (is.factor(x)) as.numeric(x) else x)
   sample_cor_before <- group

   

   for(i in 1:ncol(cor_before)){
     b <- sample_cor_before[,i]
     
     #pairwise.complete.obs
     cor_before[,i] <- cor(before, b, use = "p")
     cor_p_before[,i] <- corPvalueStudent(cor_before[,i],nSamples = nsample)
     
     #pairwise.complete.obs
     cor_after[,i] <- cor(after, b, use = "p")
     cor_p_after[,i] <- corPvalueStudent(cor_after[,i],nSamples = nsample)
   }
   

  
  
  # plot histogram of p value distribution: follow Fig. S1 in the Cell paper
  # density plot
  #windows(10,3)
  #pdf(file='cor_before_after.pdf',height=8,width=12)
  par(mfrow=c(2,3))
  # p value
  
  #number of plots to generate
  N = ncol(cor_p_after)
  
  # Insert the right number of plot output objects into the web page
  output$cor_densityPlot <- renderUI({
    
    plot_output_list <- lapply(1:N, function(i) {
      plotname <- colnames(cor_p_after)[i]
      plotlyOutput(plotname)
    })
    
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_output_list)
  })
  
  # Call renderPlot for each one. Plots are only actually generated when they
  # are visible on the web page.
  for (i in 1:N) {
    # Need local so that each item gets its own number. Without it, the value
    # of i in the renderPlot() will be the same across all instances, because
    # of when the expression is evaluated.
    local({
      my_i <- i
      plotname <- colnames(cor_p_after)[my_i]
      
      output[[plotname]] <- renderPlotly({
        
        
        #prepare dataset
        cor_p_bf <- cor_p_before[,my_i]
        cor_p_af <- cor_p_after[,my_i]
        cor_d_before <- density(cor_p_before[,my_i],na.rm = TRUE)
        cor_d_after <- density(cor_p_after[,my_i],na.rm = TRUE)
        ylim <- max(cor_d_before$y,cor_d_after$y)

        kstest <- ks.test(cor_p_bf,cor_p_af)
        ksp <- kstest$p.value

        p <- plot_ly()
        p <- add_trace(p,
                       x = ~cor_d_before$x,
                       y = ~cor_d_before$y,
                       type = 'scatter',
                       mode = 'lines',
                       name = "before",
                       line = list(color = "black")
        )
        p <- add_trace(p,
                       x = ~cor_d_after$x,
                       y = ~cor_d_after$y,
                       type = 'scatter',
                       mode = 'lines',
                       name = "after",
                       line = list(color = "red"))
        p <- p %>% layout(
          title = plotname,
          xaxis = list(title = 'p values',range = c(0,1), showgrid = F),
          yaxis = list(title = 'Density',range = c(0,ylim+0.1), showgrid = F),
          annotations = list(text = paste("p < ",ksp),  x = 0.5, y = ylim,showarrow=FALSE )
        )
        
      })
    })
  }


  
  # output$cor_densityPlot <- renderPlotly({
  #   #created tabs to put the graph
  #   N = ncol(cor_p_after)
  #   plot_list = vector("list", N)
  #   lab_list = vector("list", N)
  #   
  #   titlex = c()
  #   x = 0.15
  #   titlex = c(titlex,x)
  #   
  #   if(N == 1){
  #     titlex = x
  #   }else{
  #     for(i in 1:(N-1)){
  #       y = x + 0.4
  #       x = y
  #       titlex = c(titlex, y)
  #     }
  #   }
  #   
  #   
  #   for (k in 1:ncol(cor_p_after)) {
  #     cor_d_after <- density(cor_p_after[,k],na.rm = TRUE)
  #     cor_d_before <- density(cor_p_before[,k],na.rm = TRUE)
  #     p <- plot_ly()
  #     p <- add_trace(p,
  #                    x = ~cor_d_after$x, 
  #                    y = ~cor_d_after$y, 
  #                    type = 'scatter',
  #                    mode = 'lines',
  #                    name = "after")
  #     p <- add_trace(p,
  #                    x = ~cor_d_before$x,
  #                    y = ~cor_d_before$y,
  #                    type = 'scatter',
  #                    mode = 'lines',
  #                    name = "before"
  #                    )
  #     
  # 
  #     titley = 1.05
  #     
  #     plot_list[[k]] = p[k]
  #     lab_list[[k]] = list(x=titlex[k], y=titley, text=colnames(cor_p_after)[k], 
  #                          showarrow=F, xref='paper', yref='paper', font=list(size=18))
  # 
  #   }
  #   subplot(plot_list)%>%layout(annotations = lab_list)
    
 
    
      # p <- p %>% layout(xaxis = list(title = 'P value of correlation'),
      #                   yaxis = list(title = 'Density'))
      # p
      # plot(density(cor_p_after[,j],na.rm = TRUE),main= colnames(sample_cor_before)[j],ylim=c(0,5), xlim=c(0,1),xlab='P value of correlation',col='blue') # before
      # lines(density(cor_p_before[,j],na.rm = TRUE),col='red') # after

  # })

  
  updateProgressBar(
    session = session,
    id = "correction",
    title = "Data saved",
    value = 100
  )
  
  corRun$corRunValue <- input$runCorrection
  
  closeSweetAlert(session = session)
  sendSweetAlert(
    session = session,
    title = "DONE",
    text = "Correction Done",
    type = "success"
  )
})


output$covResultTable <- renderUI({
  if(corRun$corRunValue){
    tagList(fluidRow(column(
      12,
      downloadButton("downloadCovResult", "Download corrected table")
    )),
    tags$br(),
    fluidRow(column(
      12, 
      DT::dataTableOutput('cor_resultTable') %>% withSpinner()
    )))} else {
      helpText("Click [Run Covariates Correction] to obtain Result Table.")
    }
})

# output$cov_densityPlotUI <- renderUI({
#   if(corRun$corRunValue){
#     fluidRow(column(
#       12,
#       plotlyOutput('cor_densityPlot')%>% withSpinner()
#     ))} else {
#       helpText("Click [Run Covariates Correction] to obtain correlation density graph.")
#     }
# })
