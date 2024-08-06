
############################################################
########### Directly Source Five Native Scripts ############
############################################################
source('JUMPn_Helpers/JUMPn_functions/wgcna_mod_test_large.R')       
source('JUMPn_Helpers/JUMPn_functions/JUMPn_rich.R')      
source('JUMPn_Helpers/JUMPn_functions/PPI_Analysis.R')    
source('JUMPn_Helpers/JUMPn_functions/Network_Analysis.R')
source('JUMPn_Helpers/JUMPn_functions/JUMPn_Errors.R')    
source('JUMPn_Helpers/JUMPn_functions/ErrorModals.R')    
source('JUMPn_Helpers/JUMPn_functions/JUMPn_UIs.R') 
source('JUMPn_Helpers/JUMPn_functions/JUMPn_Modals.R')

############################################################

###################################################
###### Necessary Modules For JUMPn Execution ######
###################################################
library(shiny)                         ## Module 1
library(ggplot2)                       ## Module 2
library(DT)                            ## Module 3
library(tools)                         ## Module 4
library(tidyverse)                     ## Module 5
library(plyr)                          ## Module 6
library(writexl)                       ## Module 8
library(readxl)                        ## Module 9
library(igraph)                        ## Module 11
require(visNetwork, quietly = TRUE)    ## Module 12
library(ggnewscale)                    ## Module 13
library(WGCNA)                         ## Module 14
library(stringr)                       ## Module 15
library(annotate)                      ## Module 16
library(shinyjs)                       ## Module 17
library(shinyWidgets)
library(gplots)
####################################################

options(shiny.maxRequestSize = 30*1024^2)


########################################################
#### Define the Server Logic For Execution of JUMPn ####
########################################################
observeEvent(input$ViewExampleFile2,{
  showModal(imageModal(output))
})

observeEvent(input$ViewExampleFile3,{
  showModal(imageModal2(output))
})

observeEvent(input$FileInfo,{
  showModal(imageModal(output))
})

observeEvent(input$PPI_FileInfo,{
  showModal(imageModal2(output))
})

observeEvent(input$TransformInfo,{
  my_message <- 'By selecting "Yes", your expression data will be Log2-transformed. Log-transformation is necessary for intensity-like values (e.g., TMT reporter intensity, LFQ-intensity, iBAC protein intensity, as well as normalized spectral counts NSAF).'
  WidgetInfoModal(session,my_message)
})

observeEvent(input$FileTypeInfo,{
  my_message <- 'JUMPn currently supports .xlsx, .csv, and .txt input formats. Specifying the correct file type improves efficiency of reading the input. If you are not sure about
      features of the input file, click the "View Example File" button below.'
  WidgetInfoModal(session,my_message)
})

observeEvent(input$FileTypeInfo,{
  my_message <- 'JUMPn currently supports .xlsx, .csv, and .txt input formats. Specifying the correct file type improves efficiency of reading the input. If you are not sure about
      features of the input file, click the "View Example File" button below.'
  WidgetInfoModal(session,my_message)
})

observeEvent(input$UserFolderInfo,{
  my_message <- 'No special characters (only alphabets, numbers and underscores allowed). If an existing output folder is detected, a random number will be suffixed to the user-specified folder name to avoid overwriting.'
  WidgetInfoModal(session,my_message)
})

observeEvent(input$MinClusterInfo,{
  #my_message <- 'This parameter represents how many genes/proteins must belong to the same expression profile in order to be treated like a co-expression cluster. Otherwise, this pseudo-cluster will be merged with a neighboring cluster provided it is within the minimum similarity threshold, or it will be dropped entirely.'
  my_message <- 'Minimum number of proteins for each co-expression cluster.'
  WidgetInfoModal(session,my_message)
})

observeEvent(input$MinKMEInfo,{
  #my_message <- 'In order to identify co-expression clusters, WGCNA uses a measure of similarity to assess the strength of relationship between the expression profiles of individual genes. In JUMPn, we use Pearson Correlation as the measure of similarity. If the similarity between expression profiles exceeds this threshold, these genes will be assigned to the same cluster. This metric is applied iteratively to all genes.'
  my_message = 'Cluster membership based on Pearson correlation between individual expression trend and the eigen gene (i.e. consensus of the co-expression cluster). Lowering this value results in larger clusters. For large-sample size analysis (e.g., n > 30 samples), a small value (e.g., 0.2) is recommended.'
  WidgetInfoModal(session,my_message)
})

observeEvent(input$MinDistInfo,{
  #my_message <- 'Some expression profiles will not have enough genes to exceed the minimum cluster size specified. If this profile does however exceed the Pearson correlaion threshold with a nearby cluster, then it can be merged with that cluster.'
  my_message = 'For merging clusters after initial clustering: minimum distance (meansured by Pearson R) among eigengenes of co-expression clusters. Increasing this values results in fewer (smaller number of clusters) but larger (more proteins for each cluster on average) clusters.'
  WidgetInfoModal(session,my_message)
})

observeEvent(input$SkipPPIInfo,{
  my_message <- 'With the addition of our step-by-step analysis which allows users to separately undertake WGCNA or PPI analysis, this function will be deprecated soon.
      If you wish to undertake a different type of analysis, we suggest navigating to the Change Analysis Mode button.'
  WidgetInfoModal(session,my_message)
})

observeEvent(input$TomSizeInfo,{
  my_message = 'Modules larger than this treshold will be further split into smaller ones based on topological overlap matrix.'
  WidgetInfoModal(session,my_message)
})

observeEvent(input$TomMinInfo,{
  #my_message <- 'This parameter controls the minimum protein module size'
  my_message <- 'Modules smaller than this threshold will be discarded.'
  WidgetInfoModal(session,my_message)
})

observeEvent(input$BackgroundFileInfo,{
  my_message <- 'The gene background represents the pool of possible genes from which a set may be derived. For a pathway analysis result that better reflects your own dataset, you
      can provide your own background or merge your background with our in-house gene background.'
  WidgetInfoModal(session,my_message)
})

observeEvent(input$PPIFileInfo,{
  my_message <- 'The customized PPI file contains three tab-seperated columns in the format of: <Protein_A>, “connection”, and <Protein_B>; here <Protein_X> are presented by the official gene names of each protein.'
  WidgetInfoModal(session,my_message)
})

observeEvent(input$AnnotationFileInfo,{
  my_message <- 'As in the case of the PPI database, we also built a composite database for gene set annotation. Our composite database includes KEGG, GO, Reactome, and HALLMARK. If you wish
      to supply your own database or merge your database with ours, you may do so here. For precise examples of the annotation file format, please visit our HELP page.'
  WidgetInfoModal(session,my_message)
})

observeEvent(input$BackgroundSelectInfo,{
  my_message <- 'Here you can select the default background, a user-supplied background, or merge backgrounds'
  WidgetInfoModal(session,my_message)
})

observeEvent(input$AnnotationSelectInfo,{
  my_message <- 'Here you can select any combination of annotations'
  WidgetInfoModal(session,my_message)
})

observeEvent(input$DatabaseInfo,{
  #my_message <- 'Here ypu can choose between any of our available databases. If you wish to use your own PPI database, select the user-defined option and upload the file in the "Advanced Parameters" '
  my_message='Multiple versions of PPI databases are available for the user. Customized PPI database can be uploaded via "Advanced Parameters".'
  WidgetInfoModal(session,my_message)
})

observeEvent(input$MetaFileInfo,{
  #my_message <- 'Upload a meta file associating each sample with a specific group for large sample analysis and additional visualization.'
  my_message='A meta file specifies the sample group information. 
The meta file should contain at least two columns: column 1 must contain the sample names identical to the column names and order from the protein quantification matrix file; 
Column 2 onwards will be used for group assignment for any number of features defined by the user. 
The number of columns is flexible.'
  WidgetInfoModal(session,my_message)
})

observeEvent(input$ppi_input_file_info,{
  my_message <- 'The PPI input file must be three columns. If you already ran WGCNA, you can  use the co_exp_cluster_3columns.xlsx. Otherwise, column 1 must be group assignment for the protein, column 2 is miscellaneous infomation, and column 3 mustbe gene names as Gene Symbols'
  WidgetInfoModal(session,my_message)
})

observeEvent(input$FAQs,{
  if (input$FAQs == '1'){
    showModal(imageModal(output))
  } else if (input$FAQs == '2'){
    showModal(imageModal2(output))
  } else if (input$FAQs == '3'){
    showModal(imageModal3(output))
  } else if (input$FAQs == '4'){
    showModal(imageModal4(output))
  } else if (input$FAQs == '5'){
    showModal(imageModal5(output))
  } else if (input$FAQs == '6'){
    my_message <- 'WGCNA identifies clusters of genes/proteins that exhibit a common expression/abundance profile across replicates using several statistical parameters,
          several of which are available to the user. The minimum cluster size controls how large an identified expression/abundance profile must be to be interpreted as a cluster. 
          increasing the minimum cluster size will yield larger co-expression clusters, but small co-expression clusters may be missed.Pearson correlation reflects the degree of 
          correlatio that must be observed betweeen the profiles of individual genes/proteins for them to exist in the same cluster. We traditionally hold the Pearson Correlation at .7. 
          Decreasing the correlation will yield larger clusters, while increasing the correlation yields smaller clusters with fewer false positives. The minimum distance parameter reflects the minimum
          distance between clusters for them to be merged together. Increasing the miinimum distance will merge more clusters, and thus increase the cluster size.'
    WidgetInfoModal(session,my_message)
  } else if (input$FAQs == '7'){
    my_message <- 'Protein-Protein modularization is undertaken with two modularization steps. Initially, the protein-protein interaction network is overlaid onto the co-expression
          clusters and subsequently analyzed to identify "islands" of disconnected components. If any of these disconnected components exceed the minimum TOM-triggered module size (maximum modulle parameter for the user),
          then the second modularization step is undertaken using TOM-based modularization. In TOM-based modularization, an adjacency matrix is gennerated and analyzed to find regions of interactions that are 
          denser than the interactions throughout the network as a whole. The outcome of PPI analysis is a methodolgy which significantly enhances the identification of biological pathways and ontologies'
    WidgetInfoModal(session,my_message)
  } else if (input$FAQs == '8'){
    showModal(imageModal6(output))
  } else if (input$FAQs == '9'){
    showModal(imageModal7(output))
  }
})

hideTab('whichtab', target='WGCNAPage')
hideTab('whichtab', target='PPIPage')
the_tab <- reactive({input$whichtab})

###################################
###### Check File Extensions ######
###################################
unresolved_error_ctr <- 0

observeEvent(input$meta_file,{
  if (is.null(input$input_table)==TRUE){
    sendSweetAlert(
      session=session,
      type = 'warning',
      title='Warning:',
      btn_labels = NA,
      closeOnClickOutside = TRUE,
      text = tags$div('Please upload a WGNCA input file before providing a meta file. This way we can verify your meta file is the correct format.')
    )
    verify_input_file <- 1
  } else {
    verify_input_file <- 0
  }
  if (is.null(input$meta_file)==FALSE & verify_input_file==0){
    extension <- file_ext(input$input_table$datapath)
    input_file <- input$input_table$datapath
    if (extension == 'xlsx'){                                                                            ### If the User File is .xlsx ...
      table <- read_excel(input_file, na='NA') 
    } else if (extension == 'csv'){                                                                     ### If the User file is .csv  ...
      table <- read.csv(input_file)                                                               ### Also call read_excel function
    } else if (extension == 'txt'){                                                                     ### If the User file is ,txt
      table <- read_tsv(input_file)                                                               ### Call base r read.table function
    }                 
    col_check <- colnames(table)
    val <- Check_file_type(session,input$meta_file$datapath,c('txt','csv','xlsx'),'meta_file',col_check)
    unresolved_error_ctr <- val
  }
  #print(unresolved_error_ctr)
})

observeEvent(input$background_file,{
  if (is.null(input$background_file)==FALSE){
    print(input$background_file$datapath)
    val <- Check_file_type(session,input$background_file$datapath,c('txt','sif'),'background')
    unresolved_error_ctr <- val
  }
  #print(unresolved_error_ctr)
})
observeEvent(input$ppi_file,{
  if (is.null(input$ppi_file)==FALSE){
    print(input$ppi_file$datapath)
    val <- Check_file_type(session,input$ppi_file$datapath,c('txt','csv','xlsx','sif'),'ppi_db')
    unresolved_error_ctr <- val
  }
  #print(unresolved_error_ctr)
})
observeEvent(input$annotation_file,{
  if (is.null(input$annotation_file)==FALSE){
    print(input$annotation_file$datapath)
    val <- Check_file_type(session,input$annotation_file$datapath,c('txt','csv','xlsx'),'pathway_db')
    unresolved_error_ctr <- val
  }
  #print(unresolved_error_ctr)
})
observeEvent(input$ppi_input_file,{
  if (is.null(input$ppi_input_file)==FALSE){
    print(input$ppi_input_file$datapath)
    val <- Check_file_type(session,input$ppi_input_file$datapath,c('txt','csv','xlsx'),'ppi_input_file')
    unresolved_error_ctr <- val
  }
  #print(unresolved_error_ctr)
})
#############################################################################################
#############################################################################################
#############################################################################################

output$DynamicUI <- renderUI({
  Condensed_UI()
})

Data <- reactiveValues(
  Info = NULL
)

demo_data_ctr <<- reactive('0')

observeEvent(input$Advanced, {
  showModal(inputModal())
})

observeEvent(input$ModeofAnalysis,{
  if (input$ModeofAnalysis == "Independent PPI"){
    NewPPISelected(session)
  }
})


observeEvent(input$total_search,{

  if (is.null(input$ModeofAnalysis)==TRUE){
    ModeOfAnalysis <<- reactive("WGCNA+PPI")
  } else {
    ModeOfAnalysis <<- reactive(input$ModeofAnalysis)
  }

  KEGG <- as.data.frame(read_excel('JUMPn_Helpers/database_dependencies/ComprehensivePathwayDatabases.xlsx', sheet='KEGG'))
  Reactome <- as.data.frame(read_excel('JUMPn_Helpers/database_dependencies/ComprehensivePathwayDatabases.xlsx', sheet='Reactome'))
  GO <- as.data.frame(read_excel('JUMPn_Helpers/database_dependencies/ComprehensivePathwayDatabases.xlsx', sheet='GO'))
  HALLMARK <- as.data.frame(read_excel('JUMPn_Helpers/database_dependencies/ComprehensivePathwayDatabases.xlsx', sheet='HALLMARK'))
  
  #print(input$PathwayEnrichment)

  if (is.null(input$PathwayEnrichment)==TRUE){
    pathway_mode <- 'default'
  } else {
    pathway_mode <- input$PathwayEnrichment
  }
  
  if (pathway_mode!='User'){
    pathway_datasets <- list(KEGG, Reactome, GO, HALLMARK)
    pathway_vectorized_list <- list()
    for (db in 1:length(pathway_datasets)){
      the_db <- pathway_datasets[[db]]
      pathway_list <- strsplit(the_db[,2], '/')
      pathway_vectorized_list[[db]] <- pathway_list
    }
  } else if (pathway_mode=='User'){
    print(class(user_db))
    pathway_datasets <- list()
    pathway_datasets[[1]] <- as.data.frame(user_db)
    pathway_vectorized_list <- list()
    the_db <- as.data.frame(user_db)
    pathway_vectorized_list[[1]] <- as.list(strsplit(the_db[,2], '/'))
  }

  if (is.null(input$background_file)==TRUE){
    background_file <- read.csv("JUMPn_Helpers/database_dependencies/hs_mm_background_v3.0_20170516.txt")
  } else {
    background_file <- input$background_file
    background_file <- read.table(background_file$datapath)
  }
      #############################################################################################################################################
      #############################################################################################################################################
      ##################################################### Processing the User Input Params ######################################################
      #############################################################################################################################################
      #############################################################################################################################################
      executed_wgcna <- 0
      withProgress(message="Processing User Input File", value=0, max=1.2, {
        
        if (ModeOfAnalysis()=='WGCNA+PPI' || ModeOfAnalysis()=='Independent WGCNA'){
          initial_time <- Sys.time()                                                                        ### Get Program Start Time
          
          project_wd <- getwd()  
          project_wd
          print("line349")

          table <- variables$enrichedDataList
          idx <- grep("p.value|p-value|FDR|Log2Fold", colnames(table))
          table <- table[,-idx]
          print(table)
          print("line353")
          
          gene_names <- table[,2]                                                                           ### Collect the gene names, which will always be in column 2 
          input_data_frame <- data.frame(GN=table[,1])                                                     
          expression_data <- table[,3:ncol(table)]
          expression_title <- colnames(expression_data)
          input_data_frame <- cbind(input_data_frame, expression_data)
          
          column_name_list <- c('Gene Names')
          
          for (replicate in 1:(ncol(table)-3)){
            replicate_name <- sprintf('%d', replicate)
            column_name_list <- c(column_name_list, replicate_name)
          }
          
          colnames(input_data_frame) <- column_name_list
          
          incProgress(message = "Transforming Expression Data", .1)
          
          GN <- input_data_frame[,1]
          normalized_TMT_df <- data.frame(GN)
          input_expression_data <- as.data.frame(input_data_frame[,2:ncol(input_data_frame)])
          colnames(input_expression_data) <- colnames(table[,4:ncol(table)])
          
          transformed <- input_expression_data
          
          normalized_TMT_df <- cbind(normalized_TMT_df, transformed)
          print(normalized_TMT_df)
          #target_folder <- input$UserFolder
          target_folder <- paste0("test", str_remove_all(Sys.time(), ":"))
          dir.create(sprintf('%s/%s', project_wd, target_folder))
          
          #target_folder <- sprintf('Output_Locus/%s',target_folder)
          #dir.create(sprintf('%s/%s', project_wd, target_folder))
          write.csv(normalized_TMT_df, file=paste(target_folder, "/norm_exp_matrix.csv", sep=''), quote=T,row.names=F,sep="\t")
          #############################################################################################################################################
          #############################################################################################################################################
          #############################################################################################################################################
          #############################################################################################################################################
          
          #############################################################################################################################################
          #############################################################################################################################################
          ##################################################### Execute Co-Expression Clustering ######################################################
          #################################################################### WGCNA ##################################################################
          #############################################################################################################################################
          incProgress(message='Executing WGCNA Clustering', .1)
          min_cluster_size <- input$min_cluster_size
          min_cluster_dist <- input$min_cluster_dist
          min_kME <- input$min_kme
          network_type <- 'signed'
          scale_free_R_sq <- .8
          non_scale_free_beta <- 16
          cluster_file_list <- c()
          cluster_name_list <- c()
          enrichment_list <- list()
          execute_WGCNA <- 'Yes'
          wgcna_directory <- sprintf('%s/%s/wgcna_output',project_wd,target_folder)
          dir.create(sprintf('%s/%s/wgcna_output',project_wd,target_folder))
          pre_wgcna <- Sys.time()
          meta_file <- input$meta_file
          execute_wgcna(input_exp_matrix = sprintf('%s/norm_exp_matrix.csv',target_folder), col_offset = 1, min_cluster_size = sprintf('%.2f',min_cluster_size), min_cluster_dist = sprintf('%.2f',min_cluster_dist),  min_kME=sprintf('%.2f',min_kME), network_type=sprintf('%s',network_type), scale_free_R_sq=sprintf('%.2f',scale_free_R_sq), non_scale_free_beta=16, outputDir=sprintf('%s',wgcna_directory), meta_file)
          post_wgcna <- Sys.time()
          file_list <- list.files(path=(sprintf('%s', wgcna_directory)))
          for (item in file_list){ 
            if (substr(item, 1, 7)=='Cluster'){
              cluster_file_list <- c(cluster_file_list, item)
            }
          }
          #############################################################################################################################################
          #############################################################################################################################################
          #############################################################################################################################################
          #############################################################################################################################################
          
          #############################################################################################################################################
          #############################################################################################################################################
          ##################################################### Cluster Based Pathway Enrichment ######################################################
          #############################################################################################################################################
          #############################################################################################################################################
          
          cluster_dt_meta <- list()
          meta_ctr <- 1
          clust_ct <- 1
          cluster_specific_table_list <- list()
          dfs_for_heatmaps <- list()
          cluster_pathway_list <- list()
          cluster_specific_group_plot <- list()
          clustered_df_list <- list()
          pre_cluster_pathway_analysis <- Sys.time()
          KEGG_list <- list()
          Reactome_List <- list()
          GO_List <- list() 
          HALLMARK_List <- list()
          user_list <- list()
          
          for (file in cluster_file_list){
            search_file <- sprintf("%s/%s", wgcna_directory, file)
            cluster_data <- read.table(search_file, sep = '\t', header=TRUE)
            accession <- as.vector(unlist(cluster_data[1]))
            print(accession)
            #genes <- accession
            genes <- table[(is.element(as.vector(unlist(table[,1])),unique(accession))==TRUE),2]
            print("line454")
            print(genes)
            print("line456")
            #genes <- list(genes)
            enrichment_list <- c(enrichment_list, list(genes)) # GN list by combining all clusters; question: is this one list? Or is it N list, one for each cluster? It appears it's N list, one for each cluster!
            
            cluster_genes <- unique(toupper(unlist((genes))))
            write.table(cluster_genes, sprintf('%s/cluster%d_unique_gene_list.txt',wgcna_directory, clust_ct))
            
            pre_time <- Sys.time()
            print('pre enrich')
            if (pathway_mode!='User'){
              grouped_pathways <- pathway_enrichment2(pathway_datasets,pathway_vectorized_list, background_file,cluster_genes,5)
            } else {
              grouped_pathways <- use_defined_pathway_enrichment2(pathway_datasets,pathway_vectorized_list, background_file,cluster_genes,5)
            }
            print('post enrich')
            if (!is.null(grouped_pathways)==TRUE){
              if (!is.null(grouped_pathways[[1]])==TRUE ){
                if(nrow(grouped_pathways[[1]])>0){
                  grouped_pathway_for_export <- grouped_pathways[[1]]
                  cluster_pathway_list[[clust_ct]] <- cbind(Cluster=sprintf('Cluster %d', clust_ct),grouped_pathway_for_export)
                }
              } 
              if (!is.null(grouped_pathways[[2]])==TRUE){
                if(nrow(grouped_pathways[[2]])>0){
                  grouped_pathway_for_heatmap <- cbind(Cluster=sprintf('Cluster %d', clust_ct),grouped_pathways[[2]])
                  dfs_for_heatmaps[[clust_ct]] <- grouped_pathway_for_heatmap
                  meta_ctr <- meta_ctr + 1
                }
              }
              if (!is.null(grouped_pathways[[3]])==TRUE){
                grouped_data_list <- grouped_pathways[[3]]
                if (pathway_mode!='User'){
                  if (length(grouped_data_list)>0){KEGG_list[[clust_ct]] <- cbind(Cluster=sprintf('Cluster %d', clust_ct),Database='KEGG', grouped_data_list[[1]])}
                  if (length(grouped_data_list)>1){Reactome_List[[clust_ct]] <- cbind(Cluster=sprintf('Cluster %d', clust_ct), Database='Reactome',grouped_data_list[[2]])}
                  if (length(grouped_data_list)>2){GO_List[[clust_ct]] <- cbind(Cluster=sprintf('Cluster %d', clust_ct),Database='GO',grouped_data_list[[3]])}
                  if (length(grouped_data_list)>3){HALLMARK_List[[clust_ct]] <- cbind(Cluster=sprintf('Cluster %d', clust_ct),Database='HALLMARK',grouped_data_list[[4]])}
                } else if (pathway_mode=='User') {
                  print('made it  to the close one')
                  if (length(grouped_data_list)>0){user_list[[clust_ct]] <- cbind(Cluster=sprintf('Cluster %d', clust_ct),Database='User_db', grouped_data_list[[1]])}
                  print('error at the close one')
                }
                group_specific_data <- group_specific_enrichment(grouped_pathways[[3]],wgcna_directory,sprintf('Cluster%d', clust_ct))
                cluster_specific_table_list[[clust_ct]] <- group_specific_data[[2]]
                cluster_specific_group_plot[[clust_ct]] <- group_specific_data[[1]]
              }
            }
            table_df <- as.data.frame(table)
            df_subset <- table_df[is.element((unlist(table_df[,1])),accession)==TRUE,]
            print(df_subset)
            expression_data <-df_subset[,4:ncol(df_subset)]
            mean_vector <- (colMeans(df_subset[,4:ncol(df_subset)]))^10
            expression_labels <- colnames(df_subset[,4:ncol(df_subset)])
            expression_data <- expression_data^10
            df2 <- data.frame(labels=expression_title)
            df2 <- cbind(df2, value=as.vector(unlist(mean_vector)))
            colnames(df2) <- c('Label', 'Expression')
            df2$Label <- factor(df2$Label, levels = df2$Label)
            min <- min(expression_data)
            max <- max(expression_data)
            mean_min <- min(mean_vector) - (min(mean_vector)*.1)
            mean_max <- max(mean_vector)
            lower_range <- min - (min*.1)
            upper_range <- max + (max*.01)
            mid <- (lower_range+upper_range)/2
            cluster_title <- sprintf('Cluster %d Average Expression Profile',clust_ct)
            
            ggplot(data=df2) + 
              geom_col(aes(x=Label, y=Expression, fill=Expression)) +
              coord_cartesian(ylim=c(mean_min,mean_max)) + 
              ggtitle(cluster_title) + 
              xlab('Expression Replicates') + 
              ylab('Intensity') + 
              scale_color_continuous(low="red", high="green") + 
              theme(plot.title = element_text(size=15, face='bold', hjust=0.5),
                    axis.title.x = element_text(face='bold',size=14),
                    axis.title.y = element_text(face='bold',size=14),
                    panel.border = element_rect(colour = "black", fill=NA, size=2),
                    axis.text.x=element_text(angle=45,hjust=1,vjust=1))
            ggsave(sprintf('%s/Cluster%d_ExpressionProfile.pdf',wgcna_directory,clust_ct))
            print("line528")
            cluster_dt_meta[[clust_ct]] <- cbind(Cluster=sprintf('Cluster %d', clust_ct),df_subset)
            
            clust_ct <- clust_ct + 1
          }  
          
          cluster_title_list <- c()
          
          for (n in 1:(clust_ct-1)){
            item <- sprintf('Cluster %d', n)
            cluster_title_list <- c(cluster_title_list, item)
          }
          print("line539")
          n_clusters <- length(cluster_title_list)
          dfs_list  <- list(do.call(rbind, KEGG_list),do.call(rbind, Reactome_List),do.call(rbind, GO_List),do.call(rbind, HALLMARK_List))
          total_cluster_pathway_dataframe <- do.call(rbind,dfs_list)
          total_pathway_df <- do.call(rbind,cluster_dt_meta)
          total_cluster_pathway_df <- do.call(rbind,clustered_df_list)
          total_cluster_df <- as.data.frame(do.call(rbind,cluster_dt_meta))
          combined_cluster_pathway_table <- as.data.frame(do.call(rbind,cluster_pathway_list))
          heatmap_df <- do.call(rbind,dfs_for_heatmaps)
          pathway_heatmap(na.omit(heatmap_df),wgcna_directory,n_clusters)
          heat_file <- sprintf("%s/FDR_heatmap.png", wgcna_directory)
          #############################################################################################################################################
          #############################################################################################################################################
          #############################################################################################################################################
          #############################################################################################################################################
          write_xlsx(total_pathway_df[,1:3],sprintf('%s/co_exp_clusters_3columns.xlsx',wgcna_directory))
          executed_wgcna <- executed_wgcna + 1
          ###################################################################################################################################################
          ###################################################################################################################################################
          ################################################### Construct Co-Expression Output Objects ########################################################
          ###################################################################################################################################################
          ###################################################################################################################################################
          wgcna_files <- list.files(wgcna_directory)
          eigengene_grouplist <- c()
          eigengene_file_list <- c()
          for (file in 1:length(wgcna_files)) {
            if (substr(wgcna_files[file], 1,9)=='eigengene'){
              eigengene_file_list <- c(eigengene_file_list, wgcna_files[file])
              group <- substr(wgcna_files[file],14,nchar(wgcna_files[file])-4)
              eigengene_grouplist <- c(eigengene_grouplist, paste(group, 'Boxplot', sep = ' '))
            }
          }
          post_cluster_pathway_analysis <- Sys.time()
          if (is.null(input$meta_file)==TRUE){
            output$ExpressionDisplayOptions <- renderUI({
              selectInput('ExpressionDisplay', h5("Select the Expression Format"),
                          choices=list("Boxplot" = 'Boxplot',
                                       "Trends" = "Trends"),
                          selected = "Boxplot")
            })
          } else if (is.null(input$meta_file)==FALSE){
            output$ExpressionDisplayOptions <- renderUI({
              selectInput('ExpressionDisplay', h5("Select the Expression Format"),
                          choices=eigengene_grouplist,
                          selected = eigengene_grouplist[1])
            })
          }
          
          output$ClusterDTSelectionBox <- renderUI({
            #selectInput('ClusterDTSelect', h5("Select the Cluster You'd Like to Display"),
            selectInput('ClusterDTSelect', h5("Select the Co-expression Cluster"),
                        choices = cluster_title_list,
                        selected = 'Cluster 1')
          })
          output$DynamicClusterDTDisplay <- renderUI({
            table_title <- htmlOutput("clustertabletitle")
            alpha2 <- dataTableOutput('ClusterDT')
            cluster_content_beta <- list(table_title,alpha2)
            do.call(tagList, cluster_content_beta)
            
          })
          output$ResultsPath <- renderUI({
            path_tag <- htmlOutput('pathtag')
            do.call(tagList, list(path_tag))
          })
          the_path <- sprintf("WGCNA Data Exported To: %s \n The 3-column  PPI analysis  input file (co_exp_cluster_3columns.xlsx) may be found here.",wgcna_directory)
          path_css <- sprintf('<b style="font-size:16px; text-align:centr;">%s </b>', the_path)
          output$pathtag <- renderText(path_css)
          observeEvent(input$ClusterDTSelect, {
            print('a')
            cluster <- input$ClusterDTSelect
            #cluster_pathway_table <- total_cluster_pathway_df[total_cluster_pathway_df[,1]==cluster,]
            print('b')
            #print(combined_cluster_pathway_table[combined_cluster_pathway_table[,1]==cluster,])
            #print(nrow(combined_cluster_pathway_table[combined_cluster_pathway_table[,1]==cluster,]))
            if (nrow(combined_cluster_pathway_table[combined_cluster_pathway_table[,1]==cluster,]) >= 1){
              title_1 <- sprintf("Co-Expression %s Membership", cluster)
              css <- sprintf('<p style="font-size:30px">%s </p>', title_1)
              output$clustertabletitle <- renderText(css)
              title_d <- sprintf("Co-Expression %s Pathway Annotations Ranked", cluster)
              css2 <- sprintf('<p style="font-size:30px">%s </p>', title_d)
              output$ClusterPathwayTableTitle <- renderText(css2)
              cluster_pathway_table <- as.data.frame(combined_cluster_pathway_table[combined_cluster_pathway_table[,1]==cluster,])
              cluster_pathway_table[,4] <- formatC(as.vector(unlist(cluster_pathway_table[,4])), format = 'e',digits=2)
              cluster_pathway_table[,5] <- formatC(as.vector(unlist(cluster_pathway_table[,5])), format = 'e',digits=2)
              output$ClusterPathwayTable <- renderDataTable({
                DT::datatable(cluster_pathway_table,options = list(pageLength=10), rownames=FALSE)
              })
            } else {
              title_d <- sprintf("No Pathways Identified in %s", cluster)
              css2 <- sprintf('<p style="font-size:30px">%s </p>', title_d)
              output$ClusterPathwayTableTitle <- renderText(css2)
              empty <- data.frame('Pathway/GO', 'Overlapped Genes', 'P-Value', 'BH FDR', '# Genes in Pathway', '# Genes in Cluster/Module','# Genes Overlapped')
              output$ClusterPathwayTable <- renderDataTable({
                DT::datatable(cluster_pathway_table,options = list(pageLength=10), rownames=FALSE)
              })
            }
            print('c')
            i <- match(cluster,cluster_title_list)
            dt <- cluster_dt_meta[[i]]
            dt <- dt[,2:ncol(dt)]
            print('d')
            dt[,4:ncol(dt)] <- round(dt[,4:ncol(dt)],2)
            dt <- as.data.frame(dt)
            print(nrow(dt))
            if (nrow(dt) > 1){
              output$ClusterDT <- renderDataTable({
                DT::datatable(dt, options = list(pageLength=5))
              })
              print('e')
              observeEvent(input$ClusterDT_cell_clicked, {
                selection <- input$ClusterDT_cell_clicked    
                row_clicked <- selection$row
                if (is.null(row_clicked)==FALSE){
                  gene <- dt[row_clicked,2]
                  expression_data <- dt[row_clicked,4:ncol(table)]
                } else {
                  gene <- dt[1,2]
                  expression_data <- dt[1,4:ncol(table)]
                }
                print('f')
                mean_vector <- (colMeans(dt[,4:ncol(dt)]))
                expression_labels <- colnames(dt[,4:ncol(dt)])
                expression_data <- expression_data
                df <- data.frame(labels=expression_title)
                df <- cbind(df, value=as.vector(unlist(expression_data)))
                print('g')
                colnames(df) <- c('Label', 'Expression')
                df$Label <- factor(df$Label, levels = df$Label)
                df2 <- data.frame(labels=expression_title)
                df2 <- cbind(df2, value=as.vector(unlist(mean_vector)))
                colnames(df2) <- c('Label', 'Expression')
                df2$Label <- factor(df2$Label, levels = df2$Label)
                print('h')
                min <- min(expression_data)
                max <- max(expression_data)
                mean_min <- min(mean_vector) - (min(mean_vector)*.1)
                mean_max <- max(mean_vector)
                lower_range <- min - (min*.1)
                upper_range <- max + (max*.01)
                mid <- (lower_range+upper_range)/2
                print('i')
                title <- sprintf('%s: %s Expression Profile',cluster,gene)
                cluster_title <- sprintf('%s Average Expression Profile',cluster)
                print('j')
                fluidRow(
                  column(6,
                         output$GeneExpression <- renderPlot({
                           ggplot(data=df) + 
                             geom_col(aes(x=Label, y=Expression, fill=Expression)) +
                             coord_cartesian(ylim=c(lower_range,upper_range)) + 
                             ggtitle(title) + 
                             xlab('Expression Replicates') + 
                             ylab('Intensity') + 
                             scale_color_continuous(low="red", high="green") + 
                             theme(plot.title = element_text(size=15, face='bold', hjust=0.5),
                                   axis.title.x = element_text(face='bold',size=14),
                                   axis.title.y = element_text(face='bold',size=14),
                                   panel.border = element_rect(colour = "black", fill=NA, size=2),
                                   axis.text.x=element_text(angle=45,hjust=1, vjust=1))
                           #ggsave(sprintf('%s/Cluster%d_%s_ExpressionProfile.pdf',wgcna_directory,i,gene))
                         })),
                  print('k'),
                  column(6,
                         output$Cluster_Expression <- renderPlot({
                           ggplot(data=df2) + 
                             geom_col(aes(x=Label, y=Expression, fill=Expression)) +
                             coord_cartesian(ylim=c(mean_min,mean_max)) + 
                             ggtitle(cluster_title) + 
                             xlab('Expression Replicates') + 
                             ylab('Intensity') + 
                             scale_color_continuous(low="red", high="green") + 
                             theme(plot.title = element_text(size=15, face='bold', hjust=0.5),
                                   axis.title.x = element_text(face='bold',size=14),
                                   axis.title.y = element_text(face='bold',size=14),
                                   panel.border = element_rect(colour = "black", fill=NA, size=2),
                                   axis.text.x=element_text(angle=45,hjust=1,vjust=1))
                           #ggsave(sprintf('%s/Cluster%d_ExpressionProfile.pdf',wgcna_directory,i))
                         })),
                )
              })
            }
          })
          if (is.null(input$meta_file)==TRUE){
            observeEvent(input$ExpressionDisplay, { 
              if (input$ExpressionDisplay=='Trends'){
                output$ExpressionDisplayObject <- renderUI({
                  main_title <- htmlOutput("TrendTitle")
                  main_plot <- imageOutput("TrendPlot")
                  trend_content <- list(main_title, main_plot)
                  do.call(tagList, trend_content)
                })
                title1 <- "WGCNA Co-Expression Cluster Trends"
                css1 <- sprintf('<p style="font-size:30px">\t%s</p>', title1)
                output[['TrendTitle']] <- renderText(css1)
                trend_file <- sprintf("%s/trends_coexpression_cluster.png", wgcna_directory)
                output[['TrendPlot']] <- renderImage({
                  outfile <- tempfile(fileext='.png')
                  png(outfile, width=450, height=500)
                  dev.off()
                  list(src = trend_file,
                       contentType = 'image/png',
                       width = 450,
                       height = 500,
                       alt = "This is alternate text")
                }, deleteFile=FALSE)
              } else if (input$ExpressionDisplay=='Boxplot') {
                output$ExpressionDisplayObject <- renderUI({
                  main_title <- htmlOutput("BoxplotTitle")
                  main_plot <- imageOutput("Boxplot")
                  trend_content <- list(main_title, main_plot)
                  do.call(tagList, trend_content)
                })
                title2 <- "Boxplot Co-Expression Clusters"
                css2 <- sprintf('<p style="font-size:30px">\t%s</p>', title2)
                output[['BoxplotTitle']] <- renderText(css2)
                boxplot_file <- sprintf('%s/boxplot_coexpression_cluster.png', wgcna_directory)
                output[['Boxplot']] <- renderImage({
                  outfile <- tempfile(fileext='.png')
                  png(outfile, width=450, height=500)
                  dev.off()
                  list(src = boxplot_file,
                       contentType = 'image/png',
                       width = 450,
                       height = 500,
                       alt = "This is alternate text")
                }, deleteFile=FALSE)
              }
              
            })
          } else {
            observeEvent(input$ExpressionDisplay, { 
              num <- match(input$ExpressionDisplay,eigengene_grouplist)
              choice <- eigengene_file_list[num]
              print(choice)
              output$ExpressionDisplayObject <- renderUI({
                main_title <- htmlOutput("metaobjecttitle")
                main_plot <- imageOutput("metaobjecttplot")
                trend_content <- list(main_title, main_plot)
                do.call(tagList, trend_content)
              })
              title1 <- sprintf("Eigengene vs %s", input$ExpressionDisplay)
              css1 <- sprintf('<p style="font-size:30px">\t%s</p>', title1)
              output[['metaobjecttitle']] <- renderText(css1)
              output[['metaobjecttplot']] <- renderImage({
                outfile <- tempfile(fileext='.png')
                png(outfile, width=450, height=500)
                dev.off()
                list(src = sprintf('%s/%s', wgcna_directory,choice),
                     contentType = 'image/png',
                     width = 450,
                     height = 500,
                     alt = "This is alternate text")
              }, deleteFile=FALSE)
              
            })
          }
          #############################################################################################################################################
          #############################################################################################################################################
          ####################################################### Cluster Based Pathway Heatmap #######################################################
          #############################################################################################################################################
          #############################################################################################################################################
          output$ClusterPathwayHeatmap <- renderUI({
            alpha <- htmlOutput('HeatTitle')
            beta <- imageOutput('HeatImage')
            heat_content <- list(alpha,beta)
            do.call(tagList,heat_content)
          })
          title3 <- "Cluster Pathway Enrichment Heatmap"
          css3 <- sprintf('<p style="font-size:30px">\t%s</p>', title3)
          output[['HeatTitle']] <- renderText(css3)
          output[['HeatImage']] <- renderImage({
            outfile <- tempfile(fileext='.png')
            png(outfile, width=450, height=500)
            dev.off()
            list(src = heat_file,
                 contentType = 'image/png',
                 width = 450,
                 height = 500,
                 alt = "This is alternate text")
          }, deleteFile=FALSE)
        }
        #############################################################################################################################################
        #############################################################################################################################################
        #############################################################################################################################################
        #############################################################################################################################################
        ###################################################################################################################################################
        ###################################################################################################################################################
        ###################################################################################################################################################
        ###################################################################################################################################################
        
        
        
        #############################################################################################################################################
        #############################################################################################################################################
        ######################################################## Protein Interacton Analysis ########################################################
        #############################################################################################################################################
        #############################################################################################################################################
        print("line834")
        if (ModeOfAnalysis()=='WGCNA+PPI' || ModeOfAnalysis()=='Independent PPI'){
          print(ModeOfAnalysis())
          print("line848___")
          modular_content_list <- list()
          modular_content_names <- c()
          modular_content_list_ctr <- 1
          cluster_content_list  <- list()
          cluster_content_names <- c()
          cluster_content_list_ctr <- 1
          incProgress(message='Accessing Protein Interaction Data', .1)
          incProgress(message="Overlaying Protein Interaction Network (This May Take A Few Minutes)", .1)
          pre_ppi_analysis <- Sys.time()
          ppi_db <- input$PPIDatabases
          if (input$PPIDatabases=='User-Provided Database'){
            ppi_db <- input$ppi_file
          }
          ppi_content <- read.table('JUMPn_Helpers/database_dependencies/BioPlex_String400_Inweb150_v1.0.1.sif', header=TRUE, colClasses= c(NA, "NULL", NA))
          combined_df_list <- list()
          
          if (ModeOfAnalysis()=='WGCNA+PPI'){
            upper_combined_gene_list <- toupper(unique(unlist(enrichment_list))) # GN list by combining all clusters
            interaction_cluster <- ppi_content[((is.element(ppi_content[,1],upper_combined_gene_list)==TRUE) & (is.element(ppi_content[,2],upper_combined_gene_list)==TRUE) & (ppi_content[,1]!=ppi_content[,2])),] # subnetwork by overlaying cluster genes onto the PPI net
            for (cluster in 1:length(enrichment_list)){ # loop for each cluster: get edges within each WPC (stored in 'combined_df_list')
              combined_df_list[[cluster]] <- interaction_cluster[((is.element(interaction_cluster[,1],toupper(enrichment_list[[cluster]]))==TRUE) & (is.element(interaction_cluster[,2],toupper(enrichment_list[[cluster]]))==TRUE) ),]
              print(cluster)
            }
          } else if (ModeOfAnalysis()=='Independent PPI'){
            project_wd <- getwd() 
            #target_folder <- input$UserFolder
            if (is.element('Output_Locus', list.files(project_wd))==FALSE){
              dir.create(sprintf('%s/Output_Locus', project_wd))
            }
            print('yes')
            if (is.element(target_folder, list.files(sprintf('%s/Output_Locus',project_wd)))==TRUE){
              unlink(sprintf('%s/Output_Locus/%s', project_wd,target_folder), recursive = TRUE)
            }
            target_folder <- sprintf('Output_Locus/%s',target_folder)
            dir.create(sprintf('%s/%s', project_wd, target_folder))
            #ppi_file <- input$ppi_input_table
            ppi_file <- input$ppi_input_file
            file_type <- file_ext(ppi_file$datapath)
            if (file_type == 'xlsx'){                                                                            ### If the User File is .xlsx ...
              ppi_stuff <- read_excel(ppi_file$datapath)                                                               ### Call read_excel function on the user input file
            } else if (file_type == 'csv'){                                                                     ### If the User file is .csv  ...
              ppi_stuff <- read.csv(ppi_file$datapath)                                                               ### Also call read_excel function
            } else if (file_type == 'txt'){                                                                     ### If the User file is ,txt
              ppi_stuff <- read_tsv(ppi_file$datapath)                                                               ### Call base r read.table function
            }       
            clusters <- unique(unlist(ppi_stuff[,1]))
            for (cluster in 1:length(clusters)){
              enrichment_list <- toupper(unique(unlist(ppi_stuff[ppi_stuff[,1]==clusters[cluster],3])))
              combined_df_list[[cluster]] <-  ppi_content[((is.element(ppi_content[,1],enrichment_list)==TRUE) & (is.element(ppi_content[,2],enrichment_list)==TRUE) & (ppi_content[,1]!=ppi_content[,2])),] # subnetwork by overlaying cluster genes onto the PPI net
            }
          }
          
          print('here')
          nodular_dataframe_list <- list()
          total_cluster_ctr <- 1
          interaction_clusters <- list()
          
          modularization_folder <- 'ppi_modularization_output'
          dir.create(sprintf('%s/%s/%s',project_wd,target_folder,modularization_folder))
          
          for (cluster in combined_df_list){ # loop for each cluster: 
            print(sprintf('Cluster %d', total_cluster_ctr))
            data <- PPI_Modularization(cluster,input$TOM_triggered_module_size,total_cluster_ctr,project_wd,target_folder)
            if (is.null(data)==FALSE){
              if  (is.null(data[[1]])==FALSE){
                if (is.null(data[[2]])==FALSE){
                  nodular_dataframe_list[[total_cluster_ctr]] <- data[[1]]
                  interaction_clusters[[total_cluster_ctr]] <- data[[2]]
                  total_cluster_ctr <- total_cluster_ctr + 1
                }
              }
            }
          }
          print('there')
          #############################################################################################################################################
          #############################################################################################################################################
          #############################################################################################################################################
          #############################################################################################################################################
          
          
          post_ppi_analysis <- Sys.time()
          ###################################################################################################################################################
          ###################################################################################################################################################
          ################################################### Develop Comprehensive Network Objects #########################################################
          ###################################################################################################################################################
          ###################################################################################################################################################
          incProgress(message="Generating Total Network", .1)
          edge_frame_list <- list()
          network_data <- do.call(rbind, nodular_dataframe_list)              ##
          nodular_meta_data <- as.data.frame(network_data[,1:5])
          print("line921")
          colnames(nodular_meta_data) = c('id', 'label','group','shape','color')
          print("line923")
          final_node_frame <- nodular_meta_data[!duplicated(nodular_meta_data[,1]),]
          edge_meta_data <- do.call(rbind, interaction_clusters)
          check_list <- unlist(nodular_meta_data[,1]) 
          rownames(edge_meta_data) <- c(1:nrow((edge_meta_data)))                        ##
          colnames(edge_meta_data) <- c('from', 'to')  
          final_edge_frame <- edge_meta_data[((is.element(edge_meta_data[,'from'],check_list)==TRUE) & (is.element(edge_meta_data[,2],check_list)==TRUE) ),]
          final_edge_frame <- cbind(final_edge_frame, weights=1)
          final_edge_frame <- cbind(final_edge_frame, color='black')
          ###################################################################################################################################################
          ###################################################################################################################################################
          ###################################################################################################################################################
          ###################################################################################################################################################
          
          ###################################################################################################################################################
          ###################################################################################################################################################
          #################################################### Develop Modularized Network Objects ##########################################################
          ###################################################################################################################################################
          ###################################################################################################################################################
          incProgress(message="Creating Modularized Networks", .1)
          directory <- sprintf('%s/%s/network_models',project_wd,target_folder)
          output$PPIResultsPath <- renderUI({
            path_tag <- htmlOutput('pathtag2')
            do.call(tagList, list(path_tag))
          })
          the_path <- sprintf("PPI Data Exported To: %s",directory)
          path_css2 <- sprintf('<b style="font-size:16px; text-align:centr;">%s </b>', the_path)
          output$pathtag2 <- renderText(path_css2)
          dir.create(sprintf('%s/%s/network_models',project_wd,target_folder))
          net_file <- 'total_network.html'
          nodular_frame_list <- list() 
          modular_pathway_plot_list <- list()
          group_specific_table_list <- list()
          second_pathway_list_ctr <- 1
          modules_used <- c()
          edge_weight_list <- list()
          modular_list_ctr <- 1
          pathway_list_ctr <- 1
          pathway_df_list <- list()
          modularized_gene_list <- list()
          nod_ctr  <- 1
          cluster_modules <- list()
          cluster_module_ctr <- 1
          pre_modular_pathway_analysis <- Sys.time()
          total_modular_groups <- c()
          for(i in 1:length(nodular_dataframe_list)){
            if (!is.null(nodular_dataframe_list[[i]])==TRUE){
              data_x <- Derive_Edge_Weights(nodular_dataframe_list[[i]],interaction_clusters[[i]])
              edge_frame_list[[nod_ctr]] <- as.data.frame(data_x[[1]])
              edge_weight_list[[nod_ctr]] <- as.vector(data_x[[2]])
              nod_ctr <- nod_ctr + 1
              the_groups  <- as.vector(data_x[[4]])
              cluster_modules[[cluster_module_ctr]] <- the_groups
              cluster_module_ctr <- cluster_module_ctr + 1
              total_modular_groups <- c(total_modular_groups, the_groups)
              incProgress(message="Identifying Enriched Ontologies", .1)
              for (module in 1:length(the_groups)){
                gene_set <- unlist(data_x[[3]][data_x[[3]][,3]==the_groups[module],1])
                modularized_gene_list[[modular_list_ctr]] <- gene_set
                if (length(gene_set)>2){
                  if (pathway_mode!='User'){
                    cluster_pathway_table <- pathway_enrichment2(pathway_datasets,pathway_vectorized_list, background_file,gene_set,0)
                  } else {
                    cluster_pathway_table <- use_defined_pathway_enrichment2(pathway_datasets,pathway_vectorized_list, background_file,gene_set,0)
                  }
                  #cluster_pathway_table <- pathway_enrichment2(pathway_datasets,pathway_vectorized_list,background_file,gene_set,0)
                } else {
                  cluster_pathway_table <- NULL
                }
                if (!is.null(cluster_pathway_table[[1]])==TRUE ){
                  if(nrow(cluster_pathway_table[[1]])>0){
                    grouped_pathway_for_export <- cluster_pathway_table[[1]]
                    if (ModeOfAnalysis()=='WGCNA+PPI'){
                      cluster_pathway_list[[clust_ct]] <- cbind(Cluster=sprintf('Cluster %d', modular_list_ctr),grouped_pathway_for_export)
                    }
                  }
                } 
                if (!is.null(cluster_pathway_table[[2]])==TRUE){
                  if(nrow(cluster_pathway_table[[2]])>0){
                    grouped_pathway_for_heatmap <- cbind(Module=the_groups[module],cluster_pathway_table[[2]])
                    pathway_df_list[[pathway_list_ctr]] <- grouped_pathway_for_heatmap
                    pathway_list_ctr <- pathway_list_ctr + 1
                  }
                }
                if (!is.null(cluster_pathway_table[[3]])==TRUE){
                  group_specific_data <- group_specific_enrichment(cluster_pathway_table[[3]],directory,the_groups[module])
                  group_specific_plot_beta <- group_specific_data[[1]]
                  modular_pathway_plot_list[[second_pathway_list_ctr]] <- group_specific_plot_beta
                  group_specific_table_list[[second_pathway_list_ctr]] <- group_specific_data[[2]]
                  module_heatmap(gene_set,group_specific_data[[2]],directory,the_groups[module])
                  modules_used <- c(modules_used,the_groups[module])
                  second_pathway_list_ctr <- second_pathway_list_ctr + 1
                }
                modular_list_ctr <- modular_list_ctr + 1
              }
            }
          }
          names(modularized_gene_list) <- total_modular_groups
          post_modular_pathway_analysis <- Sys.time()
          cluster_object_list <- list()
          for (i in 1:length(nodular_dataframe_list)){
            object <- sprintf('Cluster %d', i)
            cluster_object_list[[i]] <- object
          }
          total_modules <- as.vector(unlist(cluster_modules))
          pathway_total_df <- do.call(rbind, pathway_df_list)
          modules <- unique(pathway_total_df[,1])
          lengths <- unique(pathway_total_df[,7])
          module_gene_set <- c()
          condensed_list <- c()
          size_vector <- c()
          for (module in 1:length(total_modules)){
            the_mod <- total_modules[module]
            if (the_mod %in% modules){
              pathways <- unlist(pathway_total_df[pathway_total_df[,1]==the_mod,2])
              if (length(pathways) > 10){
                pathways <- pathways[1:10]
              }
              condensed <- paste(pathways, collapse = ',')
              condensed_list <- c(condensed_list,condensed)
              module_size <- pathway_total_df[pathway_total_df[,1]==the_mod,7]
              module_gene_set <- c(module_gene_set,paste(modularized_gene_list[[the_mod]], collapse=','))
              size_vector <- c(size_vector,module_size[1])
            } else {
              condensed_list <- c(condensed_list, 'No Pathways Enriched')
              module_gene_set <- c(module_gene_set,paste(modularized_gene_list[[the_mod]], collapse=','))
              size_vector <- c(size_vector, length(modularized_gene_list[[the_mod]]))
            }
          }
          pathway_total_df <- data.frame('Module'=total_modules,"Module Proteins"=module_gene_set,'# Proteins In Module'=size_vector, 'Functional Pathways'=condensed_list)
          ###################################################################################################################################################
          ###################################################################################################################################################
          ###################################################################################################################################################
          ###################################################################################################################################################
          
          
          #####################################################################################################################################################
          #####################################################################################################################################################
          ################################################### Construct PPI Modularized Output Objects ########################################################
          #####################################################################################################################################################
          #####################################################################################################################################################
          output$DynamicClusterSelectionBox <- renderUI({
            #selectInput('ClusterSelect', h5("Select the Cluster You'd Like to Display"),
            selectInput('ClusterSelect', h5("Select the Co-expression Cluster"),
                        choices = cluster_object_list,
                        selected = 'Cluster 1')
          })
          output$DynamicClusterDisplay <- renderUI({
            alpha <- htmlOutput('ClusterPlotTitle')
            beta <- visNetworkOutput('ClusterPlot')
            gamma <- uiOutput('ClusterBorder')
            fred <- uiOutput('LayoutSelectionPanel')
            delta <- dataTableOutput('ClusterTable')
            epsilon <- dataTableOutput('ClusterCoefficient')
            cluster_content <- list(alpha,beta,gamma,fred,delta,epsilon)
            do.call(tagList, cluster_content)
          })
          output$DynamicModularDisplay <- renderUI({
            alpha3 <- htmlOutput('ModularPlotTitle')
            beta3 <- visNetworkOutput('ModularPlot')
            gamma3 <- uiOutput('ModularBorder')
            delta3 <- dataTableOutput('ModularTable')
            epsilon3 <- dataTableOutput('ModularCoefficient')
            modular_content <- list(alpha3,beta3,gamma3,delta3,epsilon3)
            do.call(tagList, modular_content)
          })
          output$LayoutSelectionPanel <- renderUI({
            selectInput('LayoutSelection', h5("Select a Network Layout Format"),
                        choices=list("Fruchterman Reingold" = 'layout_with_fr',
                                     "Circle" = "layout_in_circle"),
                        selected = "Fruchterman Reingold")
          })
          
          
          observeEvent(input$ClusterSelect, {
            cluster <- input$ClusterSelect
            i <- match(cluster,cluster_object_list)
            network_content <- nodular_dataframe_list[[i]]
            num_genes <- nrow(network_content)
            group_vector1 <- unique(network_content[,'group'])
            num_modules <- length(group_vector1)
            title_1 <- sprintf("Co-Expression Cluster %d: %d Modules Found, %d Constituent Genes", i, num_modules, num_genes)
            css <- sprintf('<p style="font-size:30px">%s </p>', title_1)
            output$ClusterPlotTitle <- renderText(css)
            layout <- reactive({ 
              input$LayoutSelection
            })
            observeEvent(layout(), {
              output$ClusterPlot <- renderVisNetwork({
                if (layout() == 'layout_with_fr'){
                  visNetwork(network_content, edge_frame_list[[i]], height='1400px', width='100%') %>% 
                    visOptions(highlightNearest = TRUE, manipulation=TRUE,selectedBy = "group") %>%
                    visIgraphLayout(layout='layout_with_fr', weights=edge_weight_list[[i]]) 
                } else if (layout() == 'layout_in_circle') { 
                  visNetwork(network_content, edge_frame_list[[i]], height='1400px', width='100%') %>% 
                    visOptions(highlightNearest = TRUE, manipulation=TRUE,selectedBy = "group") %>%
                    visIgraphLayout(layout='layout_in_circle') 
                }
              })
            })
            output$ClusterBorder <- renderUI({
              addCssClass(class = 'networkplotborder1', selector = 'ClusterPlot')
            })
            visNetwork(network_content, edge_frame_list[[i]], main=title_1 , height='2000px', width='100%') %>% 
              visOptions(highlightNearest = TRUE, manipulation=TRUE) %>%
              visIgraphLayout(layout='layout_with_fr', weights=edge_weight_list[[i]]) %>%
              visSave(file=sprintf("%s/cluster_%d_network.html",directory,i), selfcontained = FALSE)
            
            output$DynamicModuleSelectionBox <- renderUI({
              selectInput('ModuleSelect', h5("Select the Module You'd Like to Display"),
                          choices = cluster_modules[[i]])
            })
            output$EMAP <- renderUI({
              alpha <- htmlOutput('ModulePathwayTitle')
              beta <- uiOutput('ModulePathwayPlot')
              the_best <- list(alpha,beta)
              do.call(tagList,the_best)
            })
            output$ModulePathwaySelection <- renderUI({
              selectInput('ModulePathwaySelect', h5("Select the Pathway Annotation Style"),
                          choices = list('Barplot'='Barplot',
                                         'Circle Plot' = 'Circle Plot',
                                         'Heatmap'='Heatmap',
                                         'Table' = 'Table'))
            })
            modulepathwayplot <- reactive({ 
              input$ModulePathwaySelect
            })
            observeEvent(input$ModuleSelect, {
              module <- input$ModuleSelect
              j <- match(module,cluster_modules[[i]])
              k <- match(module,modules_used)
              gene_frame <- final_node_frame[(final_node_frame[,'group']==module),]
              gene_vector <- as.vector(gene_frame[,'id'])
              edges <- unique(final_edge_frame[((is.element(final_edge_frame[,'from'], gene_vector)==TRUE) & (is.element(final_edge_frame[,'to'], gene_vector)==TRUE)),])
              edge_vector <- as.vector(unique(unlist(list(edges[,1],edges[,2]))))
              gene_frame <- subset(gene_frame,id %in% edge_vector  )
              edge_weights <- rep(1, nrow(edges))
              module_title <- sprintf("%s: %d Constituent Genes", module, length(gene_vector))
              css <- sprintf('<p style="font-size:30px">%s </p>', module_title)
              output$ModularPlotTitle <- renderText(css)
              output$ModularPlot <- renderVisNetwork({
                visNetwork(gene_frame, edges, height='1400px', width='100%') %>% 
                  visOptions(highlightNearest = TRUE, manipulation=TRUE) %>%
                  visIgraphLayout(layout='layout_with_fr', weights=edge_weights) 
              })
              
              observeEvent(modulepathwayplot(),{
                if (is.null(group_specific_table_list[[k]])==FALSE){
                  if (modulepathwayplot()=='Barplot') {
                    output$ModulePathwayPlot <- renderUI({
                      alpha <- plotOutput('ModulePathwayFigure')
                      the_best <- list(alpha)
                      do.call(tagList,the_best)
                    })
                    this_title <- sprintf('%s Pathway Barplot', module)
                    css <- sprintf('<p style="font-size:25px">%s </p>', this_title)
                    output$ModulePathwayTitle <- renderText(css)
                    output$ModulePathwayFigure <- renderPlot({
                      modular_pathway_plot_list[[k]]
                    },height=400,width=900)
                  } else if (modulepathwayplot()=='Table') {
                    output$ModulePathwayPlot <- renderUI({
                      alpha <- dataTableOutput('ModulePathwayTable')
                      the_best <- list(alpha)
                      do.call(tagList,the_best)
                    })
                    this_title <- sprintf('%s Pathway Table', module)
                    css <- sprintf('<p style="font-size:25px">%s </p>', this_title)
                    new_table <- group_specific_table_list[[k]]
                    new_table[,4] <- formatC(10^-(as.vector(unlist(new_table[,4]))), format='e',  digits=2)
                    output$ModulePathwayTitle <- renderText(css)
                    output$ModulePathwayTable <- renderDataTable({
                      DT::datatable(as.data.frame(new_table, options = list(pageLength=12)))
                    })
                  } else if (modulepathwayplot()=='Circle Plot') {
                    output$ModulePathwayPlot <- renderUI({
                      alpha <- imageOutput('ModulePathwayCircle')
                      the_best <- list(alpha)
                      do.call(tagList,the_best)
                    })
                    this_title <- sprintf('%s Pathway Circle Plot', module)
                    css <- sprintf('<p style="font-size:25px">%s </p>', this_title)
                    output$ModulePathwayTitle <- renderText(css)
                    output$ModulePathwayCircle <- renderPlot({
                      circular.barplot(group_specific_table_list[[k]])
                    },height = 550,width = 700)
                  } else if (modulepathwayplot()=='Heatmap') {
                    output$ModulePathwayPlot <- renderUI({
                      alpha <- imageOutput('ModuleHeatmap')
                      the_best <- list(alpha)
                      do.call(tagList,the_best)
                    })
                    this_title <- sprintf('%s Pathway Heatmap', module)
                    css <- sprintf('<p style="font-size:25px">%s </p>', this_title)
                    output$ModulePathwayTitle <- renderText(css)
                    output$ModuleHeatmap <- renderImage({
                      outfile <- tempfile(fileext='.png')
                      png(outfile, width=450, height=500)
                      dev.off()
                      list(src = sprintf("%s/%s_heatmap.png",directory,module),
                           contentType = 'image/png',
                           width = 450,
                           height = 500,
                           alt = "This is alternate text")
                    }, deleteFile=FALSE)
                    
                  }
                } else if (is.null(group_specific_table_list[[k]])==TRUE){
                  output$ModulePathwayPlot <- renderUI({
                    alpha <- imageOutput('NullImage')
                    the_best <- list(alpha)
                    do.call(tagList,the_best)
                  })
                  output$NullImage <- renderImage({
                    outfile <- tempfile(fileext='.png')
                    png(outfile, width=450, height=500)
                    dev.off()
                    list(src = sprintf("%s/JUMPn_Helpers/Null_Objects/null_pathway.png",project_wd),
                         contentType = 'image/png',
                         width = 450,
                         height = 500,
                         alt = "This is alternate text")
                  }, deleteFile=FALSE)
                }
              })
            })
            if (is.null(input$meta_file)==TRUE){
              if (ModeOfAnalysis()!='Independent PPI'){
                output$TheClusterExpressionDisplay <- renderUI({
                  selectInput('ClusterExpressionDisplay',
                              h5("Select the Expression Format"),
                              choices=list("Boxplot" = 'Boxplot',
                                           "Trends" = "Trends",
                                           "Pathway Barplot"="Pathway Barplot",
                                           "Pathway Circle Plot"="Pathway Circle Plot"),
                              selected = "Boxplot",
                              width = '50%')
                })
              }
            } else {
              #choice_vector <- eigengene_grouplist
              choice_vector <<- c("Pathway Barplot",eigengene_grouplist, "Pathway Circle Plot")
              output$TheClusterExpressionDisplay <- renderUI({
                selectInput('ClusterExpressionDisplay',
                            h5("Select the Expression Format"),
                            choices=choice_vector,
                            selected = choice_vector[1],
                            width = '50%')
              })
            }
            observeEvent(input$ClusterExpressionDisplay, {
              if (ModeOfAnalysis()=='WGCNA+PPI' || ModeOfAnalysis()=='Independent WGCNA'){
                if (is.null(input$meta_file)==TRUE){
                  if (input$ClusterExpressionDisplay=='Trends'){
                    output$ExpressionFormat <- renderUI({
                      new_main_title <- htmlOutput("NewTrendTitle")
                      new_main_plot <- imageOutput("NewTrendPlot")
                      new_trend_content <- list(new_main_title, new_main_plot)
                      do.call(tagList, new_trend_content)
                    })
                    title2 <- sprintf("Cluster %s Co-Expression Trend", i)
                    css2 <- sprintf('<p style="font-size:30px">\t%s</p>', title2)
                    output[['NewTrendTitle']] <- renderText(css2)
                    new_trend_file <- sprintf('%s/trend_coexpression_cluster_%d.png', wgcna_directory,i)
                    output$NewTrendPlot <- renderImage({
                      outfile <- tempfile(fileext='.png')
                      png(outfile, width=400, height=500)
                      dev.off()
                      list(src = new_trend_file,
                           contentType = 'image/png',
                           width = 400,
                           height = 500,
                           alt = "This is alternate text")
                    }, deleteFile=FALSE)
                  } else if (input$ClusterExpressionDisplay=='Boxplot') {
                    output$ExpressionFormat <- renderUI({
                      new_main_title <- htmlOutput("NewBoxplotTitle")
                      new_main_plot <- imageOutput("NewBoxplot")
                      new_trend_content <- list(new_main_title, new_main_plot)
                      do.call(tagList, new_trend_content)
                    })
                    title2 <- sprintf("Cluster %s Boxplot Trend", i)
                    css2 <- sprintf('<p style="font-size:30px">\t%s</p>', title2)
                    output[['NewBoxplotTitle']] <- renderText(css2)
                    new_boxplot_file <- sprintf('%s/small_boxplot_coexpression_cluster_%d.png', wgcna_directory,i)
                    output[['NewBoxplot']] <- renderImage({
                      outfile <- tempfile(fileext='.png')
                      png(outfile, width=400, height=500)
                      dev.off()
                      list(src = new_boxplot_file,
                           contentType = 'image/png',
                           width = 400,
                           height = 500,
                           alt = "This is alternate text")
                    }, deleteFile=FALSE)
                  } else if (input$ClusterExpressionDisplay=='Pathway Barplot') {
                    output$ExpressionFormat <- renderUI({
                      main_title <- htmlOutput("PathwayTitle")
                      main_plot <- plotOutput("PathwayBarPlot")
                      new_trend_content <- list(main_title, main_plot)
                      do.call(tagList, new_trend_content)
                    })
                    title2 <- sprintf("Cluster %d Annotations", i)
                    css2 <- sprintf('<p style="font-size:30px">\t%s</p>', title2)
                    output[['PathwayTitle']] <- renderText(css2)
                    output[['PathwayBarPlot']] <- renderPlot({
                      cluster_specific_group_plot[[i]]
                    }, height = 500, width=900)
                    
                  } else if (input$ClusterExpressionDisplay=='Pathway Circle Plot') {
                    output$ExpressionFormat <- renderUI({
                      main_title <- htmlOutput("PathwayTitle")
                      main_plot <- plotOutput("PathwayCirclePlot")
                      new_trend_content <- list(main_title, main_plot)
                      do.call(tagList, new_trend_content)
                    })
                    title2 <- sprintf("Cluster %d Annotations", i)
                    css2 <- sprintf('<p style="font-size:30px">\t%s</p>', title2)
                    output[['PathwayTitle']] <- renderText(css2)
                    output[['PathwayCirclePlot']] <- renderPlot({
                      circular.barplot(cluster_specific_table_list[[i]])
                    }, height = 500, width=800)
                    
                  }
                } else {
                  num <- match(input$ClusterExpressionDisplay,choice_vector)
                  choice <- eigengene_file_list[num]
                  print(choice)
                  if ( (num != length(choice_vector)) & (num !=  1) ){
                    output$ExpressionFormat <- renderUI({
                      main_title <- htmlOutput(sprintf("metaobjecttitle%d",num))
                      main_plot <- imageOutput(sprintf("metaobjecttplot%d", num))
                      trend_content <- list(main_title, main_plot)
                      do.call(tagList, trend_content)
                    })
                    title1 <- sprintf("Eigengene vs %s", input$ClusterExpressionDisplay)
                    css1 <- sprintf('<p style="font-size:30px">\t%s</p>', title1)
                    output[[sprintf("metaobjecttitle%d",num)]] <- renderText(css1)
                    output[[sprintf("metaobjecttplot%d", num)]] <- renderImage({
                      outfile <- tempfile(fileext='.png')
                      png(outfile, width=450, height=500)
                      dev.off()
                      list(src = sprintf('%s/%s', wgcna_directory,choice),
                           contentType = 'image/png',
                           width = 450,
                           height = 500,
                           alt = "This is alternate text")
                    }, deleteFile=FALSE)
                  } else if (num == 1){
                    output$ExpressionFormat <- renderUI({
                      main_title <- htmlOutput("metaobjecttitle1")
                      main_plot <- imageOutput("metaobjecttplot1")
                      trend_content <- list(main_title, main_plot)
                      do.call(tagList, trend_content)
                    })
                    title1 <- sprintf("Cluster %d Annotations", i)
                    css2 <- sprintf('<p style="font-size:30px">\t%s</p>', title1)
                    output[['metaobjecttitle1']] <- renderText(css2)
                    output[['metaobjecttplot1']] <- renderPlot({
                      cluster_specific_group_plot[[i]]
                    }, height = 500, width=900)
                  } else if (num == length(choice_vector)){
                    output$ExpressionFormat <- renderUI({
                      main_title <- htmlOutput("metaobjecttitle0")
                      main_plot <- imageOutput("metaobjecttplot0")
                      trend_content <- list(main_title, main_plot)
                      do.call(tagList, trend_content)
                    })
                    title1 <- sprintf("Cluster %d Pathway Circle Plot", i)
                    css3 <- sprintf('<p style="font-size:30px">\t%s</p>', title1)
                    output[['metaobjecttitle0']] <- renderText(css3)
                    output[['metaobjecttplot0']] <- renderPlot({
                      circular.barplot(cluster_specific_table_list[[i]])
                    }, height = 500, width=800)
                  }
                  
                  
                }
              } else if (ModeOfAnalysis()=='Independent PPI'){
              }
            }) 
          })
          incProgress(message="Computing Network Modularity", .1)
          incProgress(message="Creating Modularized Networks", .1)
          group_list <- list()
          lapply(1:length(nodular_dataframe_list), function(i){
            if (!is.null(nodular_dataframe_list[[i]])==TRUE){
              network_content <- nodular_dataframe_list[[i]]
              num_genes <- nrow(network_content)
              group_vector1 <- unique(network_content[,'group'])
              num_modules <- length(group_vector1)
              title_1 <- sprintf("Co-Expression Cluster %d: %d Modules Found, %d Constituent Genes", i, num_modules, num_genes)
              visNetwork(network_content, edge_frame_list[[i]], main=title_1 , height='2000px', width='100%') %>% 
                visOptions(highlightNearest = TRUE, manipulation=TRUE) %>%
                visIgraphLayout(layout='layout_with_fr', weights=edge_weight_list[[i]]) %>%
                visSave(file=sprintf("%s/cluster_%d_network.html",directory,i), selfcontained = FALSE)
            }
          })
        }
        #####################################################################################################################################################
        #####################################################################################################################################################
        #####################################################################################################################################################
        #####################################################################################################################################################
        log_file <- sprintf('
      JUMPn Version 1.0.0
      Data Exported to: %s/%s
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ~~ User Input Parameters ~~
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Min Cluster Size: %f 
      Min Pearson Correlation: %f 
      Min Distance: %f 
      Min Module Size: %f 
      PPI Database: %s
      ', project_wd, target_folder, 
                            input$min_cluster_size,input$min_kme, input$min_cluster_dist,
                            input$TOM_triggered_module_size, input$PPIDatabases)
        
        write(log_file, sprintf('%s/%s/JUMPn.log',project_wd,target_folder))
        file_path <- sprintf('%s/%s/ComprehensiveSummaryTables.xlsx',project_wd,target_folder)
        #observeEvent(input$downloadData, {
        if (ModeOfAnalysis()=='WGCNA+PPI'){
          publication_data <- write_xlsx(list("Cluster Genes"=as.data.frame(total_pathway_df),
                                              "Cluster Enriched Pathways"=as.data.frame(total_cluster_pathway_dataframe),
                                              "Module Enriched Pathways"=as.data.frame(pathway_total_df)),
                                         path=file_path)
          write_xlsx(list("Cluster Genes"=as.data.frame(total_pathway_df),
                          "Cluster Enriched Pathways"=as.data.frame(total_cluster_pathway_dataframe),
                          "Module Enriched Pathways"=as.data.frame(pathway_total_df)),
                     path=file_path)
        } else if (ModeOfAnalysis()=='Independent WGCNA'){
          publication_data <- write_xlsx(list("Cluster Genes"=as.data.frame(total_pathway_df),"Cluster Enriched Pathways"=as.data.frame(total_cluster_pathway_dataframe)),
                                         path=file_path)
          write_xlsx(list("Cluster Genes"=as.data.frame(total_pathway_df),"Cluster Enriched Pathways"=as.data.frame(total_cluster_pathway_dataframe)),
                     path=file_path)
        } else if (ModeOfAnalysis()=='Independent PPI'){
          print(as.data.frame(pathway_total_df))
          publication_data <- write_xlsx(as.data.frame(pathway_total_df),path=file_path)
          write_xlsx(as.data.frame(pathway_total_df),path=file_path)
        }
        #successful_download(session)
        #})
        print(list.files(sprintf('%s/%s', project_wd,target_folder)))
        output_directory <<- sprintf('%s/%s', project_wd,target_folder)
        
        
        incProgress(message="Program Finished! Go To Network Analysis Page to View Results", .1)
      })
      total_time <- Sys.time()
      # print(initial_time)
      # print(pre_wgcna)
      # print(post_wgcna)
      # print(post_cluster_pathway_analysis)
      # print(post_cluster_pathway_analysis)
      # print(pre_ppi_analysis)
      # print(post_ppi_analysis)
      # print(pre_modular_pathway_analysis)
      # print(post_modular_pathway_analysis)
      # print(total_time)
      # 
      # Total_Run_Time <- total_time-initial_time
      # WGCNA_time <- post_wgcna-pre_wgcna
      # cluster_pathway_time <- post_cluster_pathway_analysis-pre_cluster_pathway_analysis
      # ppi_time <- post_ppi_analysis - pre_ppi_analysis
      # modular_pathway_time <- post_modular_pathway_analysis - pre_modular_pathway_analysis
      # total_pathway_time <- cluster_pathway_time + modular_pathway_time
      # print(Total_Run_Time)
      # print(WGCNA_time)
      # print(ppi_time)
      # print(cluster_pathway_time)
      # print(modular_pathway_time)
      # print(total_pathway_time)
      print(gc())
      #the_locale <- sprintf('%s/%s', project_wd,input$UserFolder)
      programFinishedModal(session,output_directory)
      observeEvent(input$SubmitQuickPPI, {
        ModeOfAnalysis()=='WGCNA+PPI'
      })
      
      if (ModeOfAnalysis()=='WGCNA+PPI'){
        showTab('whichtab', target='WGCNAPage')
        showTab('whichtab', target='PPIPage')
        output$WGCNATitlePage<- renderText({
          'Results Page 1: WGCNA Output'
        })
        output$PPITitlePage<- renderText({
          'Results Page 2: PPI Output'
        })
      } else if (ModeOfAnalysis()=='Independent WGCNA'){
        showTab('whichtab', target='WGCNAPage')
        output$WGCNATitlePage<- renderText({
          'Results Page: WGCNA Output'
        })
      } else if (ModeOfAnalysis()=='Independent PPI'){ 
        showTab('whichtab', target='PPIPage')
        output$PPITitlePage<- renderText({
          'Results Page: PPI Output'
          
        })
      }
      
      observeEvent(input$TakeMeThere,{
        if (ModeOfAnalysis()=='WGCNA+PPI'){
          updateTabsetPanel(session, 'whichtab',selected='WGCNAPage')
          closeSweetAlert()
          ArrivedToResultsModal('WGCNA+PPI',session)
        }else if (ModeOfAnalysis()=='Independent WGCNA'){
          updateTabsetPanel(session, 'whichtab',selected='WGCNAPage')
          closeSweetAlert()
          ArrivedToResultsModal('Independent WGCNA',session)
          observeEvent(input$SubmitQuickPPI,{
            removeModal()
            showModal(QuickPPIModal())
          })
        } else if (ModeOfAnalysis()=='Independent PPI'){
          updateTabsetPanel(session, 'whichtab',selected='PPIPage')
          closeSweetAlert()
          ArrivedToResultsModal('Independent PPI',session)
        }
      })
      observeEvent(input$ViewResults,{
        closeSweetAlert()
      })
      berify <- reactive(input$TotalData)
      print(input$TotalData)
      observeEvent(berify(),{
        if (berify()==TRUE){
          print(berify())
          output$downloadSelectedData1 <- downloadHandler(
            filename = function() {
              paste('JUMPn_Output_Data', ".zip", sep = "")
            },
            content = function(file) {
              files2zip <- dir(sprintf('%s/%s',project_wd, target_folder), full.names = TRUE)
              print(files2zip)
              zip(zipfile=file, files=files2zip)    
            }
          )
          successful_download(session)
        }
        
      })
    })

