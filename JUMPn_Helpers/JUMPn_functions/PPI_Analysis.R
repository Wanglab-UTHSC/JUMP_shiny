##### PPI Tools #########
#########################

###################################################################
#### Derive Edge Weights Based on Module/Cluster Relationships ####
###################################################################
PPI_Modularization <- function(Cluster_Edge_Frame, min_TOM_triggered_module_size,cluster_index,root_dir,export_loc){
  if (min_TOM_triggered_module_size < 5){
    min_TOM_triggered_module_size <- 4
  }
  colnames(Cluster_Edge_Frame) <- c('from', 'to')
  sub_cluster_list <- unique(unlist(list(Cluster_Edge_Frame[,1],Cluster_Edge_Frame[,2]))) # genes within PPI subnet of this cluster
  igraph_object <- graph_from_edgelist(as.matrix(Cluster_Edge_Frame))
  components <- decompose.graph(igraph_object)
  color_list <- c('red','yellow', 'orange', 'grey', 'pink', 'teal', 'violet','cadetblue','chartreuse','aquamarine', 'orange', 'grey','steelblue1', 'pink', 'teal', 'plum1')
  novel_cluster_df_list <- list()
  novel_edge_df_list <- list()
  novel_cluster_df_list_ctr <- 1
  module_ctr <- 1
  if (length(components)!=0 & length(sub_cluster_list > 0)){
    for (item in 1:length(components)){
      #print('round')
      nodes <- V(components[[item]])
      #print('1')
      edges2 <- Cluster_Edge_Frame[((is.element(Cluster_Edge_Frame[,1],names(nodes))==TRUE) & (is.element(Cluster_Edge_Frame[,2],names(nodes))==TRUE)),]
      #print('2')
      edges2 <- rbind(edges2,data.frame(from=edges2[,2], to=edges2[,1]))
      #print('3')
      if (length(nodes) >= min_TOM_triggered_module_size){
       # print('4')
        adjacency_matrix <- matrix(data=0, nrow=length(unique(unlist(edges2))),
                                   ncol=length(unique(unlist(edges2))),
                                   dimnames=list(sort(unique(unlist(edges2))),sort(unique(unlist(edges2)))))
        diag(adjacency_matrix) <- 1
        adjacency_matrix[as.matrix(edges2)] <- 1
        TOM = TOMsimilarity(adjacency_matrix);
        dissTOM = 1-TOM
        geneTree = hclust(as.dist(dissTOM), method = "average");
        dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                    deepSplit = 2, pamRespectsDendro = FALSE,
                                    minClusterSize = 4);
        #print('5')
        dynamicColors = labels2colors(dynamicMods)
        tmp=table(dynamicColors)
        tmp=tmp[order(tmp,decreasing=T)]
        m.colors=names(tmp)
        ctr <- 1
        #print('6')
        module_dataframe_list <- list()
        #print('7')
        for (i in 1:length(tmp)) {
          group=m.colors[i]
          gene_list <- as.vector(rownames(adjacency_matrix)[dynamicColors==group])
          if (!is.null(gene_list)==TRUE){
            #print('8')
            module_edge <- edges2[is.element(edges2[,'from'],gene_list)==TRUE & is.element(edges2[,'to'],gene_list)==TRUE,]
            #print('9')
            igraph_object <- graph_from_edgelist(as.matrix(module_edge))
            #print('10')
            components2 <- decompose.graph(igraph_object)
            #print('11')
            #print(length(components2))
            if (length(components2)==1){
              module_length <- length(gene_list)
              module_name <- sprintf("Module_%s",module_ctr)
              groupname <- sprintf('Cluster%d_Module%d', cluster_index, module_ctr)
              novel_module_dataframe <- data.frame(id=gene_list, label=gene_list, group=groupname, shape='hexagon', color=color_list[ctr])
              module_dataframe_list[[module_ctr]] <- novel_module_dataframe
              #write.table(rownames(adjacency)[dynamicColors==group],file=paste(root_dir, '/', export_loc, '/', 'ppi_modularization_output/', "Second_level_module_cluster_",cluster_index,"_", 'module_', module_ctr, group,".txt",sep=""),quote=F,col.name=F,row.name=F,sep="\t")
              ctr <- ctr + 1
              module_ctr <- module_ctr + 1
            } else if (length(components2)>1){
                for (graph_comp in 1:length(components2)){
                  nodes <- names(V(components2[[graph_comp]]))
                  #print(nodes)
                  module_length <- length(nodes)
                  module_name <- sprintf("Module_%s",module_ctr)
                  groupname <- sprintf('Cluster%d_Module%d', cluster_index, module_ctr)
                  novel_module_dataframe <- data.frame(id=nodes, label=nodes, group=groupname, shape='hexagon', color=color_list[module_ctr])
                  module_dataframe_list[[module_ctr]] <- novel_module_dataframe
                  #write.table(rownames(adjacency)[dynamicColors==group],file=paste(root_dir, '/', export_loc, '/', 'ppi_modularization_output/', "Second_level_module_cluster_",cluster_index,"_", 'module_', module_ctr, group,".txt",sep=""),quote=F,col.name=F,row.name=F,sep="\t")
                  ctr <- ctr + 1
                  module_ctr <- module_ctr + 1
                }
            }
          }
        }
        combined_module_dataframe <- do.call(rbind, module_dataframe_list)
        novel_cluster_df_list[[novel_cluster_df_list_ctr]] <- combined_module_dataframe
        novel_edge_df_list[[novel_cluster_df_list_ctr]] <- edges2
        novel_cluster_df_list_ctr <- novel_cluster_df_list_ctr + 1
        
        plotCorD = dissTOM^(1);
        diag(plotCorD) = NA;
        TOM_file <- sprintf('%s/ppi_modularization_output/TOM_matrix_%d.png', export_loc, cluster_index)
        png(file=TOM_file,width=400,height=400)
        TOMplot(plotCorD, geneTree, dynamicColors, main = "")
        dev.off()
      } else {
        #print('problem is here')
        module_name <- sprintf("Module_%s",module_ctr)
        #print('a')
        groupname <- sprintf('Cluster%d_Module%d', cluster_index, module_ctr)
        #print('b')
        new_nodes <- unique(unlist(edges2))
        #print('c')
        if (length(new_nodes >=  1)){
         # print('d')
          node_content_dataframe <- data.frame(id = new_nodes, label = new_nodes, group=groupname, shape='hexagon', color=color_list[module_ctr])
          #print('e')
          novel_cluster_df_list[[novel_cluster_df_list_ctr]] <- node_content_dataframe
          #print('f')
          novel_edge_df_list[[novel_cluster_df_list_ctr]] <- edges2
          #print('g')
          novel_cluster_df_list_ctr <- novel_cluster_df_list_ctr + 1
          #print('h')
          module_ctr <- module_ctr + 1
        }
      }
    }
  }
  #print('i')
  cluster_node_df <- do.call(rbind, novel_cluster_df_list)
  #print('j')
  cluster_edge_df <- do.call(rbind, novel_edge_df_list)
  #print('k')
  if (!is.null(cluster_node_df)==TRUE){
    return(list(cluster_node_df,cluster_edge_df))
  }
}
