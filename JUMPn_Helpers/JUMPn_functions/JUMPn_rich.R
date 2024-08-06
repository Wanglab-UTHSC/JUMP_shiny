
pathway_enrichment2 <- function(pathway_file,pathway_sep_list, background_genes,query_genes,set_limit){
  ##################################################################
  ################### Define Helper Functions ######################
  ##################################################################
  #### Get Length of Overlap #######################################                                  
  overlapped_length <- function(pathway_set,query_set){            #
    length(intersect(pathway_set,query_set))                       #
  }                                                                #
  #### Length of Overlap Between Nested Lists ######################                 
  overlapped_list <- function(pathway_set,query_set,fungshun){     #
    mapply(fungshun,pathway_set,query_set)                         #
  }                                                                #
  #### Get Overlapped Gene Names ###################################                              
  overlapped_genes <- function(pathway_set,query_set){             #
    intersect(pathway_set,query_set)                               #
  }                                                                #
  #### Get Genes Not Found in Background ###########################                      
  genes_not_in_background <- function(query_set,background){       ######
    setdiff(background,query_set)
  }                                                                ######
  ##################################################################
  Vectorized_Fisher_Function <- function(component_df){
    fisher_mat <- matrix(as.numeric(component_df[c(1:4)]),ncol=2,byrow = TRUE)
    fisher <- fisher.test(as.table(fisher_mat),alternative = 'greater')
    return(fisher$p.value)
  }
  ##########################################
  ###### Define Available Databases ########
  ##########################################
  set <- c('KEGG','Reactome','GO','HALLMARK')
  ##########################################
  ##########################################
  
  ###########################################
  #### Search Each Database Individually ####
  ###########################################
  group_pathways <- list()
  used_vector <- c()
  group_pathways_for_export <-  list()
  real_set_ctr <- 1
  for (database_source in 1:length(set)) {
    #### Grab the Database and Vectorize the Genes in Each Pathway ####
    pathway_database <- pathway_file[[database_source]]
    pathway_list <- pathway_sep_list[[database_source]]
    #### Relate the Pathway Name to the Vectorized Gene Set ####
    check_df <- cbind(pathway_database[,1],pathway_list)
    rownames(check_df) <- c(1:nrow(check_df))
    #### Call Overlap Length Function, Get Overlap Between Query Set and Every Pathway (This is Fast)
    overlap_length_vector <- sapply(pathway_list,overlapped_length,query_set=toupper(query_genes) )
    pathway_df1 <- as.data.frame(cbind(pathway_list,overlap_length_vector))
    rownames(pathway_df1) <- c(1:nrow(pathway_df1))
    #### Filter Out Pathways With Less than 3 overlapped genes
    pathway_df <- pathway_df1[pathway_df1[,2]>set_limit,]
    filtered_pathways <- pathway_df[,1]
    filtered_pathway_names <- row.names(pathway_df)
    if (length(filtered_pathways) >= 1){
      now_list <- as.vector(unlist(check_df[filtered_pathway_names,1]))
      if (is.null(now_list)==FALSE){
        overlapped_gene_list <- sapply(filtered_pathways,overlapped_genes,query_set=toupper(query_genes) )
        if (length(overlapped_gene_list) >= 1){ 
          background_not_in_path <- sapply(filtered_pathways,genes_not_in_background,background=unlist(background_genes,use.names=FALSE))
          if (typeof(background_not_in_path)=='character'){
            background_not_in_path <- split(background_not_in_path,seq(ncol(background_not_in_path)))
          }
          background_not_in_overlap <- sapply(overlapped_gene_list,genes_not_in_background,background=unlist(background_genes,use.names = FALSE))
          if (typeof(background_not_in_overlap)=='character'){
            background_not_in_overlap <- split(background_not_in_overlap,seq(ncol(background_not_in_overlap)))
          }
          pre_fisher_list <- list(overlapped_gene_list,filtered_pathways,unname(background_not_in_path),unname(background_not_in_overlap))
          #### Executing Fisher Test in Vectorized Fashion #####
          #### Compile Fisher Test Components Into a Dataframe By Row ####
          a <- sapply(pre_fisher_list[[2]],overlapped_length, toupper(query_genes))
          b <- overlapped_list(pre_fisher_list[[2]],pre_fisher_list[[4]],overlapped_length)
          c <- sapply(pre_fisher_list[[4]],overlapped_length,toupper(query_genes))
          d <- overlapped_list(pre_fisher_list[[3]],pre_fisher_list[[4]],overlapped_length)
          fisher_df <- data.frame(a,b,c,d)
          #### Call the Fisher Test Function on the Data Frame ####
          p_val <- apply(fisher_df,1, Vectorized_Fisher_Function)
          ##### Construct the  Dataframe ########
          if (typeof(overlapped_gene_list)=='character'){
            if (is.null(ncol(overlapped_gene_list))==FALSE){
              overlapped_gene_list <- split(overlapped_gene_list,seq(ncol(overlapped_gene_list)))
            } else if (is.null(ncol(overlapped_gene_list))==TRUE){
              overlapped_gene_list <- split(overlapped_gene_list,seq(length(overlapped_gene_list)))            
            }
          }
          overlapped_gene_fix <- unlist(lapply(overlapped_gene_list, paste, collapse='/'))
          overlapped_lengths <- lengths(overlapped_gene_list)
          pathway_length <- lengths(filtered_pathways)
          total_df <- data.frame('Pathway/GO'=unname(now_list), 'Overlapped Genes'=unname(overlapped_gene_fix), 'P-Value'=unname(p_val))
          total_df <- unique(cbind(total_df, 'BH FDR' = p.adjust(p_val,'BH'), '# Genes in Pathway'=pathway_length,'# Genes in Cluster/Module' = length(query_genes), '# Genes Overlapped'=unname(unlist(overlapped_lengths))))[,1:7]
          group_pathways[[real_set_ctr]] <- total_df[order(total_df[,4]),]
          group_pathways_for_export[[real_set_ctr]] <- total_df[order(total_df[,4]),]
          used_vector <- c(used_vector, set[database_source])
          real_set_ctr <- real_set_ctr + 1
        }
      }
    }
  }
  combined_df <- do.call(rbind,group_pathways)
  if (is.null(combined_df)==FALSE) {
    combined_total_df <- combined_df[order(combined_df[,4]),]
    combined_filtered_df <-  combined_total_df[combined_total_df[,4]<.15,]
    #if(!is.null(combined_filtered_df) &&  nrow(combined_filtered_df)>0){ 
    #  if (nrow(combined_filtered_df>5)){
    #    combined_filtered_df <- combined_filtered_df[1:6,]
    #  } 
    #}
    names(group_pathways_for_export) <- used_vector
    combined_total_df <- combined_total_df[combined_total_df[,4]<.2,]
    return(list(combined_total_df,combined_filtered_df,group_pathways_for_export))
  }
}

use_defined_pathway_enrichment2 <- function(pathway_file,pathway_sep_list, background_genes,query_genes,set_limit){
  ##################################################################
  ################### Define Helper Functions ######################
  ##################################################################
  #### Get Length of Overlap #######################################                                  
  overlapped_length <- function(pathway_set,query_set){            #
    length(intersect(pathway_set,query_set))                       #
  }                                                                #
  #### Length of Overlap Between Nested Lists ######################                 
  overlapped_list <- function(pathway_set,query_set,fungshun){     #
    mapply(fungshun,pathway_set,query_set)                         #
  }                                                                #
  #### Get Overlapped Gene Names ###################################                              
  overlapped_genes <- function(pathway_set,query_set){             #
    intersect(pathway_set,query_set)                               #
  }                                                                #
  #### Get Genes Not Found in Background ###########################                      
  genes_not_in_background <- function(query_set,background){       ######
    setdiff(background,query_set)
  }                                                                ######
  ##################################################################
  Vectorized_Fisher_Function <- function(component_df){
    fisher_mat <- matrix(as.numeric(component_df[c(1:4)]),ncol=2,byrow = TRUE)
    fisher <- fisher.test(as.table(fisher_mat),alternative = 'greater')
    return(fisher$p.value)
  }
  ##########################################
  ###### Define Available Databases ########
  ##########################################
  set <- c('KEGG','Reactome','GO','HALLMARK')
  ##########################################
  ##########################################
  
  ###########################################
  #### Search Each Database Individually ####
  ###########################################
  group_pathways <- list()
  used_vector <- c()
  group_pathways_for_export <-  list()
  real_set_ctr <- 1
    pathway_database <- pathway_file[[1]]
    pathway_list <- pathway_sep_list[[1]]
    #### Relate the Pathway Name to the Vectorized Gene Set ####
    check_df <- cbind(pathway_database[,1],pathway_list)
    rownames(check_df) <- c(1:nrow(check_df))
    #### Call Overlap Length Function, Get Overlap Between Query Set and Every Pathway (This is Fast)
    overlap_length_vector <- sapply(pathway_list,overlapped_length,query_set=toupper(query_genes) )
    pathway_df1 <- as.data.frame(cbind(pathway_list,overlap_length_vector))
    rownames(pathway_df1) <- c(1:nrow(pathway_df1))
    #### Filter Out Pathways With Less than 3 overlapped genes
    pathway_df <- pathway_df1[pathway_df1[,2]>set_limit,]
    filtered_pathways <- pathway_df[,1]
    filtered_pathway_names <- row.names(pathway_df)
    if (length(filtered_pathways) >= 1){
      now_list <- as.vector(unlist(check_df[filtered_pathway_names,1]))
      if (is.null(now_list)==FALSE){
        overlapped_gene_list <- sapply(filtered_pathways,overlapped_genes,query_set=toupper(query_genes) )
        if (length(overlapped_gene_list) >= 1){ 
          background_not_in_path <- sapply(filtered_pathways,genes_not_in_background,background=unlist(background_genes,use.names=FALSE))
          if (typeof(background_not_in_path)=='character'){
            background_not_in_path <- split(background_not_in_path,seq(ncol(background_not_in_path)))
          }
          background_not_in_overlap <- sapply(overlapped_gene_list,genes_not_in_background,background=unlist(background_genes,use.names = FALSE))
          if (typeof(background_not_in_overlap)=='character'){
            background_not_in_overlap <- split(background_not_in_overlap,seq(ncol(background_not_in_overlap)))
          }
          pre_fisher_list <- list(overlapped_gene_list,filtered_pathways,unname(background_not_in_path),unname(background_not_in_overlap))
          #### Executing Fisher Test in Vectorized Fashion #####
          #### Compile Fisher Test Components Into a Dataframe By Row ####
          a <- sapply(pre_fisher_list[[2]],overlapped_length, toupper(query_genes))
          b <- overlapped_list(pre_fisher_list[[2]],pre_fisher_list[[4]],overlapped_length)
          c <- sapply(pre_fisher_list[[4]],overlapped_length,toupper(query_genes))
          d <- overlapped_list(pre_fisher_list[[3]],pre_fisher_list[[4]],overlapped_length)
          fisher_df <- data.frame(a,b,c,d)
          #### Call the Fisher Test Function on the Data Frame ####
          p_val <- apply(fisher_df,1, Vectorized_Fisher_Function)
          ##### Construct the  Dataframe ########
          if (typeof(overlapped_gene_list)=='character'){
            if (is.null(ncol(overlapped_gene_list))==FALSE){
              overlapped_gene_list <- split(overlapped_gene_list,seq(ncol(overlapped_gene_list)))
            } else if (is.null(ncol(overlapped_gene_list))==TRUE){
              overlapped_gene_list <- split(overlapped_gene_list,seq(length(overlapped_gene_list)))            
            }
          }
          overlapped_gene_fix <- unlist(lapply(overlapped_gene_list, paste, collapse='/'))
          overlapped_lengths <- lengths(overlapped_gene_list)
          pathway_length <- lengths(filtered_pathways)
          total_df <- data.frame('Pathway/GO'=unname(now_list), 'Overlapped Genes'=unname(overlapped_gene_fix), 'P-Value'=unname(p_val))
          total_df <- unique(cbind(total_df, 'BH FDR' = p.adjust(p_val,'BH'), '# Genes in Pathway'=pathway_length,'# Genes in Cluster/Module' = length(query_genes), '# Genes Overlapped'=unname(unlist(overlapped_lengths))))[,1:7]
          group_pathways[[real_set_ctr]] <- total_df[order(total_df[,4]),]
          group_pathways_for_export[[real_set_ctr]] <- total_df[order(total_df[,4]),]
          used_vector <- c(used_vector, 'User Defined DB')
          real_set_ctr <- real_set_ctr + 1
        }
      }
    }
  combined_df <- do.call(rbind,group_pathways)
  if (is.null(combined_df)==FALSE) {
    combined_total_df <- combined_df[order(combined_df[,4]),]
    combined_filtered_df <-  combined_total_df[combined_total_df[,4]<.15,]
    #if(!is.null(combined_filtered_df) &&  nrow(combined_filtered_df)>0){ 
    #  if (nrow(combined_filtered_df>5)){
    #    combined_filtered_df <- combined_filtered_df[1:6,]
    #  } 
    #}
    names(group_pathways_for_export) <- used_vector
    combined_total_df <- combined_total_df[combined_total_df[,4]<.2,]
    return(list(combined_total_df,combined_filtered_df,group_pathways_for_export))
  }
}

pathway_heatmap <- function(pathway_dataframe,directory_path,number_of_clusters){
  unique_group_vector <- unique(pathway_dataframe[,1])
  unique_pathway_vector <- unique(pathway_dataframe[,2])
  FDR_matrix <- matrix(data=1e-2, nrow=length(unique_pathway_vector),
                       ncol=length(unique_group_vector),
                       dimnames=list(unique_pathway_vector,unique_group_vector))
  for (row in 1:nrow(pathway_dataframe)){
    FDR_matrix[pathway_dataframe[row,2],pathway_dataframe[row,1]] <- pathway_dataframe[row,4]
  }
  FDR_matrix=-log10(FDR_matrix)
  min_val <- min(FDR_matrix)
  pathways <- c()
  for (col in 1:ncol(FDR_matrix)){
    pathways <- c(pathways,names(FDR_matrix[order(FDR_matrix[,col],decreasing = TRUE),col])[1:5])
  }
  #print(pathways)
  FDR_matrix <- as.data.frame(FDR_matrix[unique(pathways),1:ncol(FDR_matrix)])
  n_row <- nrow(FDR_matrix)
  column_dif <- as.integer(number_of_clusters - ncol(FDR_matrix))
  #print(number_of_clusters)
  #print(column_dif)
  #print(FDR_matrix)
  #print('a')
  if (!is.null(column_dif)==TRUE){
    if (length(column_dif)==0){
      column_dif <- 0
    }
    #print(column_dif)
    if (column_dif > 0){
      #print('a.2')
      for (i in 1:column_dif){
        #print('a.3')
        title <- sprintf("Cluster %d", (ncol(FDR_matrix)+1))
        zero_vector <- rep(2,n_row)
        FDR_matrix <- cbind(FDR_matrix, zero_vector)
      }
    }
  }
  #print(FDR_matrix)
  #print('b')
  title_vector <- c()
  print(ncol(FDR_matrix))
  for (i in 1:ncol(FDR_matrix)){
    title_vector <- c(title_vector, sprintf('Cluster %d', i))
  }
  #print('c')
  colnames(FDR_matrix) <- title_vector
  #print('d')
  #print(FDR_matrix)
  #print(ncol(FDR_matrix))
  if (ncol(FDR_matrix) > 1){
    my_palette <- colorRampPalette(c( "white","blue"))(n = 299)
    heat_file <- sprintf("%s/FDR_heatmap.png", directory_path)
    par(mar=c(1,1,1,1))
    png(file=heat_file,width=1000,height=800)
    heatmap.2(as.matrix(FDR_matrix),scale='n' ,
              Colv = FALSE,
              density.info='n',
              trace="none",         # turns off trace lines inside the heat map
              margins = c(12,35),     # widens margins around plot
              col=my_palette,       # use on color palette defined earlier ,
              dendrogram ="r",
              key.xlab = '-log10(p-value)',
              colsep=0:ncol(FDR_matrix),rowsep=0:nrow(FDR_matrix),sepwidth=c(0.001,0.001),sepcolor="grey",keysize =1,srtCol=45,cexCol = 1.4,cexRow = 1.4
    )
    dev.off()
  } else if (ncol(FDR_matrix) == 1){
    my_palette <- colorRampPalette(c( "white","blue"))(n = 299)
    heat_file <- sprintf("%s/FDR_heatmap.png", directory_path)
    par(mar=c(1,1,1,1))
    png(file=heat_file,width=1000,height=1000)
    #print('here')
    #print(as.matrix(cbind(FDR_matrix,FDR_matrix[,1])))
    heatmap.2(as.matrix(cbind(FDR_matrix,FDR_matrix[,1])),
              scale='n' ,
              Colv = FALSE,
              density.info='n',
              trace="none",         # turns off trace lines inside the heat map
              margins = c(12,28),     # widens margins around plot
              col=my_palette,       # use on color palette defined earlier ,
              dendrogram ="r",
              rowsep=0:nrow(FDR_matrix),sepwidth=c(0.001,0.001),sepcolor="grey",keysize =1,srtCol=45,cexCol = 2,cexRow = 2
    )
    dev.off()
  }
  #print('e')
}
module_heatmap <- function(module_gene_set, module_pathway_dataframe,directory_path,module_title){
  pathways <- strtrim(unlist(module_pathway_dataframe[,2]),25)
  pathway_p <- unlist(module_pathway_dataframe[,4])
  pathway_genes <- lapply(module_pathway_dataframe[,3], strsplit, '/')
  genes <- module_gene_set 
  FDR_matrix <- matrix(data=1e-2, nrow=length(genes),
                       ncol=length(pathways),
                       dimnames=list(genes,pathways))
  for (pathway in 1:length(pathway_genes)) {
    pathway_set <- unlist(pathway_genes[[pathway]]) 
    FDR_matrix[pathway_set,pathways[pathway]] <- pathway_p[pathway]
  }
  if (ncol(FDR_matrix) > 1){
    my_palette <- colorRampPalette(c( "white","blue"))(n = 299)
    heat_file <- sprintf("%s/%s_heatmap.png", directory_path,module_title)
    par(mar=c(1,1,1,1))
    png(file=heat_file,width=800,height=1000)
    heatmap.2(FDR_matrix,scale='n' ,
              Colv = FALSE,
              density.info='n',
              trace="none",         # turns off trace lines inside the heat map
              margins = c(28,12),     # widens margins around plot
              col=my_palette,       # use on color palette defined earlier ,
              dendrogram ="r",
              colsep=0:ncol(FDR_matrix),rowsep=0:nrow(FDR_matrix),sepwidth=c(0.001,0.001),sepcolor="grey",keysize =1,srtCol=45,cexCol = 2,cexRow = 2
    )
    dev.off()
  }
}


group_specific_enrichment <- function(grouped_pathway_list,directory,module){
  grouped_pathway_list <- as.list(grouped_pathway_list)
  group_df_list <- list()
  for (group  in  1:length(grouped_pathway_list)){
    database  <- names(grouped_pathway_list)[group]
    if (nrow(grouped_pathway_list[[group]])>3){
      pathway_names  <- grouped_pathway_list[[group]][1:3,1]
      p_vals  <- -log10(grouped_pathway_list[[group]][1:3,3])
      overlapped_genes <- grouped_pathway_list[[group]][1:3,2]
    } else {
      pathway_names  <- grouped_pathway_list[[group]][,1]
      p_vals  <- -log10(grouped_pathway_list[[group]][,3])
      overlapped_genes <- grouped_pathway_list[[group]][,2]
    }
    group_df_list[[group]] <- data.frame(Pathway=pathway_names,OverlappedGenes=overlapped_genes,PVal=p_vals)
    group_df_list[[group]] <- cbind(Groups=database, group_df_list[[group]])
  }
  grouped_top_df <- as.data.frame(do.call(rbind, group_df_list))
  plot_y <- ggplot(data=grouped_top_df) +
    geom_col(aes(x=Pathway,y=PVal,fill=Groups))  + 
    coord_flip() +
    ggtitle('Enriched Pathways and Ontologies') +
    ylab(expression(-log[10](p-value))) +
    scale_x_discrete(label = function(x) strtrim(x, 45)) +
    theme(plot.title = element_text(face='bold'),
          axis.text.x = element_text(angle=45,hjust=1,vjust=1),
          axis.text.y = element_text(size=10,face='bold'))

  return(list(plot_y, grouped_top_df))
}

circular.barplot<-function(object_df){
  if (is.null(object_df)==FALSE) {
    if (nrow(object_df)>1) {
      par(mar=c(1,1,1,1))
      object_df[,2] <- strtrim(object_df[,2],40)
      df<-data.frame(values=sort(object_df[,4]), labels=object_df[order(object_df[,4]),2],group=object_df[order(object_df[,4]),1])
      plot(NA,xlim=c(-2,1.5),ylim=c(-1.5,1.5),axes=F, xlab=NA, ylab=NA, asp=1,frame.plot = TRUE)
      box("figure", lwd=2)
      grid()
      t<-sapply(df$values,function(x).5*pi-seq(0, x/5,length=1000))
      max_val <- max(df$values)
      convert <- df$values/max_val
      t <- sapply(convert,function(x).5*pi-seq(0, x*1.5*pi,length=1000))
      x<-sapply(1:nrow(df),function(x)(.3+x/nrow(df))*cos(t[,x]))
      y<-sapply(1:nrow(df),function(x)(.3+x/nrow(df))*sin(t[,x]))
      
      gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
      }
      cols <- gg_color_hue(4)
      
      for(i in 1:nrow(df)){
        if (df[i,3]=='GO'){
          color <- cols[1]
        } else if (df[i,3]=='HALLMARK'){
          color <- cols[2]
        } else if (df[i,3]=='KEGG'){
          color <- cols[3]
        } else if (df[i,3]=='Reactome'){
          color <- cols[4]
        }
        lines(x=x[,i],y=y[,i],col=color,lwd=12,lend=1)
        text(x[1,i],y[1,i],paste(df$labels[i]," - ",round(df$values[i],2),sep=""),
             pos=2,cex=.75)
      }
      dbs <- c('GO','HALLMARK','KEGG','Reactome')
      dbs <- dbs[1:length(unique(df$group))]
      print(unique(df$group))
      print(dbs)
      legend(-1.9,-.5,legend = dbs,col=unique(cols),pch=rep(15,length(unique(df$group))))
    }
  }
}



modular_pathway_enrichment2 <- function(pathway_file,background_genes,query_genes,set_limit){
  ##################################################################
  ################### Define Helper Functions ######################
  ##################################################################
  #### Get Length of Overlap #######################################                                  
  overlapped_length <- function(pathway_set,query_set){            #
    length(intersect(pathway_set,query_set))                       #
  }                                                                #
  #### Length of Overlap Between Nested Lists ######################                 
  overlapped_list <- function(pathway_set,query_set,fungshun){     #
    mapply(fungshun,pathway_set,query_set)                         #
  }                                                                #
  #### Get Overlapped Gene Names ###################################                              
  overlapped_genes <- function(pathway_set,query_set){             #
    intersect(pathway_set,query_set)                               #
  }                                                                #
  #### Get Genes Not Found in Background ###########################                      
  genes_not_in_background <- function(query_set,background){       ######
    background[is.element(background,query_set)==FALSE]            #
  }                                                                ######
  ##################################################################
  ##########################################
  ###### Define Available Databases ########
  ##########################################
  set <- c('KEGG','Reactome','GO','HALLMARK')
  ##########################################
  ##########################################
  
  ###########################################
  #### Search Each Database Individually ####
  ###########################################
  group_pathways <- list()
  used_vector <- c()
  group_pathways_for_export <-  list()
  real_set_ctr <- 1
  for (database_source in 1:length(set)) {
    #### Grab the Database and Vectorize the Genes in Each Pathway ####
    pathway_database <- as.data.frame(read_xlsx(pathway_file, sheet = set[database_source]))
    pathway_list <- strsplit(pathway_database[,2], '/')
    
    #### Relate the Pathway Name to its Vectorized Gene Set ####

    #### Call Overlap Length Function, Get Overlap Between Query Set and Every Pathway (This is Fast)
    overlap_length_vector <- overlapped_list(pathway_list,toupper(query_genes),overlapped_length)
    pathway_df1 <- cbind(pathway_list,overlap_length_vector)
    
    #### Filter Out Pathways With Less than 3 overlapped genes
    pathway_df <- as.data.frame(pathway_df1[pathway_df1[,2]>set_limit,])
    filtered_pathways <- pathway_df[,1]
    if (length(filtered_pathways) >= 1){
      now_list <- as.vector(unlist(check_df[is.element(check_df[,2],filtered_pathways)==TRUE,1]))
      if (is.null(now_list)==FALSE){
        overlapped_gene_list <- sapply(filtered_pathways,overlapped_genes,query_set=toupper(query_genes) )
        if (length(overlapped_gene_list) >= 1){ 
          background_not_in_path <- sapply(filtered_pathways,genes_not_in_background,background=unlist(background_genes,use.names=FALSE))
          #background_not_in_path <- rapply(background_not_in_path,na.omit,how='replace')
          background_not_in_overlap <- sapply(overlapped_gene_list,genes_not_in_background,background=unlist(background_genes,use.names = FALSE))
          if (typeof(background_not_in_overlap)=='character'){
            background_not_in_overlap <- split(background_not_in_overlap,seq(ncol(background_not_in_overlap)))
          }
          #background_not_in_overlap <- rapply(background_not_in_overlap,na.omit,how='replace')
          pre_fisher_list <- list(overlapped_gene_list,filtered_pathways,unname(background_not_in_path),unname(background_not_in_overlap))
          
          #### Executing Fisher Test in Vectorized Fashion #####
          Vectorized_Fisher_Function <- function(component_df){
            fisher_mat <- matrix(as.numeric(component_df[c(1:4)]),ncol=2,byrow = TRUE)
            fisher <- fisher.test(as.table(fisher_mat),alternative = 'greater')
            return(fisher$p.value)
          }
          
          #### Compile Fisher Test Components Into a Dataframe By Row ####
          a <- sapply(pre_fisher_list[[2]],overlapped_length, toupper(query_genes))
          b <- overlapped_list(pre_fisher_list[[2]],pre_fisher_list[[4]],overlapped_length)
          c <- sapply(pre_fisher_list[[4]],overlapped_length,toupper(query_genes))
          d <- overlapped_list(pre_fisher_list[[3]],pre_fisher_list[[4]],overlapped_length)
          fisher_df <- data.frame(a,b,c,d)
          #### Call the Fisher Test Function on the Data Frame ####
          p_val <- apply(fisher_df,1, Vectorized_Fisher_Function)
          ##### Construct the  Dataframe ########
          overlapped_gene_fix <- unlist(lapply(overlapped_gene_list, paste, collapse='/'))
          total_df <- data.frame('Pathway/GO'=now_list, 'Overlapped Genes'=overlapped_gene_fix, 'P-Value'=p_val)
          total_df <- cbind(total_df,p.adjust(p_val,'BH'))
          total_df <- total_df[order(total_df[,4]),]
          group_pathways[[real_set_ctr]] <- total_df
          group_pathways_for_export[[real_set_ctr]] <- subset(total_df, total_df[,4]<.4)
          used_vector <- c(used_vector, set[database_source])
          real_set_ctr <- real_set_ctr + 1
        }
      }
    }
  }
  combined_df <- do.call(rbind,group_pathways)
  if (is.null(combined_df)==FALSE) {
    combined_total_df <- combined_df[order(combined_df[,4]),]
    combined_filtered_df <-  combined_total_df[combined_total_df[,4]<.3,]
    if(!is.null(combined_filtered_df) &&  nrow(combined_filtered_df)>0){ 
      if (nrow(combined_filtered_df>5)){
        combined_filtered_df <- combined_filtered_df[1:6,]
      } 
    }
    names(group_pathways_for_export) <- used_vector
    combined_total_df <- combined_total_df[combined_total_df[,4]<.2,]
    return(list(combined_total_df,combined_filtered_df,group_pathways_for_export))
  }
}

ppi_based_pathway_enrichment <- function(pathway_database,background_genes,ppi_database,query_genes) {
  overlapped_length <- function(pathway_set,query_set){
    length(intersect(pathway_set,query_set))
  }
  pathway_list <- strsplit(pathway_database[,2], '/')
  overlap_length_vector <- sapply(pathway_list,overlapped_length,query_set=query_genes )
  #pathway_df1 <- cbind(pathway_list,overlap_length_vector)
  #pathway_df <- pathway_df1[pathway_df1[,2]>3,]
  #filtered_pathways <- pathway_df[,1]
  pathway_interactions <- list()
  for (pathways in 1:length(pathway_list)){
    pathway <- toupper(unlist(pathway_list[[pathways]]))
    interactions <- ppi_database[((is.element(ppi_database[,1],unlist(toupper(pathway)))==TRUE) & (is.element(ppi_database[,2],unlist(toupper(pathway)))==TRUE) & (ppi_database[,1]!=ppi_database[,2])),] # subnetwork by overlaying cluster genes onto the PPI net
    pathway_interactions[[pathways]] <- list(as.data.frame(unique(interactions)))
  }
  names(pathway_interactions) <- unlist(pathway_database[,1])
  capture.output(summary(pathway_interactions), 'PathwayEdgeDatabase.txt')
  
}
