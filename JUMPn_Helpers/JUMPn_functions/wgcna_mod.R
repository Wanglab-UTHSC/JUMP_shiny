# inputs:
outputDir='.'
#input_exp_matrix='norm_exp_matrix.txt'
#col_offset=3

#min_cluster_size=5
#min_cluster_dist=0.1
#min_kME=0.7

#network_type='signed'
#scale_free_R_sq=0.8
#non_scale_free_beta=16

#-----------------------------------------------------------------------------
# To parse R parameters
# Credit to Nurcan Tuncbag from MIT
#args=(commandArgs(TRUE))
#for (e in commandArgs(T)) {
 # ta = strsplit(e,"=",fixed=TRUE)
  #var = ta[[1]][1]
  #if(! is.na(ta[[1]][2])) {
   # temp = ta[[1]][2]
    #var = substr(ta[[1]][1],2,nchar(ta[[1]][1]))
    #if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "I") {
     # temp = as.integer(temp)
    #}
    #if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "N") {
     # temp = as.numeric(temp)
    #}
    #if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "V") {
     # temp = strsplit(temp,',')[[1]]
    #}
    #assign(var,temp)
    #cat("assigned ",var," the value of |",temp,"|\n")
  #} else {
   # var_fields = strsplit(var,'-')[[1]]
    #var = var_fields[length(var_fields)]
    #assign(var,TRUE)
    #cat("assigned ",var," the value of TRUE\n")
  #}
#}

#-----------------------------------------------------------------------------
execute_wgcna <- function(input_exp_matrix,network_type,outputDir,col_offset,min_cluster_dist,min_cluster_size,min_kME,scale_free_R_sq,non_scale_free_beta,metaFile){
  # make sure WGCNA is installed: https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/index.html#cranInstall
  
  # string to number
  outputDir = outputDir
  col_offset=as.numeric(col_offset)
  
  min_cluster_dist=as.numeric(min_cluster_dist)
  min_cluster_size=as.numeric(min_cluster_size)
  min_kME=as.numeric(min_kME)
  
  scale_free_R_sq=as.numeric(scale_free_R_sq)
  non_scale_free_beta=as.numeric(non_scale_free_beta)
  
  # prepare matrix for WGCNA
  whl=read.csv(input_exp_matrix,head=T,quote = "", fill=TRUE)
  tb=whl[,(col_offset+1):ncol(whl)]
  
  plex=ncol(tb)
  m=t(tb)
  datExpr=m
  
  # Load the WGCNA package
  library(WGCNA)
  
  # select beta
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  #powers = c(seq(from = 20, to=40, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = network_type)
  
  # print figure of beta selection
  pdf(file=paste(outputDir,'/sftPwr_select.pdf',sep=''),width=10,height=10)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.80,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
  
  # set beta so that R^2>scale_free_R_sq (default: 0.8); but if this beta > max_beta_scale_free (default: 20), set beta=non_scale_free_beta (default: 16)
  if (sum(-sign(sft$fitIndices[,3])*sft$fitIndices[,2]>scale_free_R_sq)==0) { # not scale free for beta < 20
    #if (sum(-sign(sft$fitIndices[,3])sft$fitIndices[,2]>scale_free_R_sq)==0) { # not scale free for beta < 20
    softPower = non_scale_free_beta
  } else {
    softPowerRow = which(sft$fitIndices[,2]>scale_free_R_sq)[1]
    softPower = sft$fitIndices[softPowerRow,1]
  }
  
  # adjacency matrix
  softPower
  adjacency = adjacency(datExpr, power = softPower,type = network_type);
  #adjacency = lapply(adjacency, as.numeric)
  
  # Turn adjacency into topological overlap
  TOM = TOMsimilarity(adjacency);
  dissTOM = 1-TOM
  
  # Call the hierarchical clustering function
  geneTree = hclust(as.dist(dissTOM), method = "average");
  
  # We like large modules, so we set the minimum module size relatively high:
  minModuleSize = min_cluster_size;
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize);
  table(dynamicMods)
  
  # Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  
  # Plot the dendrogram and colors underneath
  pdf(file=paste(outputDir,'/dendrogram_raw.pdf',sep=''),width=10,height=10)
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  dev.off()
  
  # normalize signal to print trend (like what STEM do)
  nm=tb[,1:plex]
  for (i in 1:plex) {
    nm[,i]=tb[,i]-tb[,1]
  }
  summary(nm)
  head(nm)
  dim(nm)
  
  tmp=table(dynamicColors);length(tmp)
  tmp=tmp[order(tmp,decreasing=T)]
  #tmp=tmp[-which(names(tmp)=='grey')] # remove 'grey' module
  
  # added on 1/23/19
  # if too many (>25) modules generated, only print the top 25
  if (length(tmp)>25) tmp=tmp[1:25]
  # added end
  
  m.colors=tmp
  png(file=paste(outputDir,'/boxplot_raw.png',sep=''),width=10000,height=10000)
  box_col=ceiling(sqrt(length(tmp)))
  box_row=ceiling(length(tmp)/box_col)
  par(las=2,mfrow=c(box_row,box_col),mar = c(1, 1, 1, 1),cex.lab=1.2)
  for (i in 1:length(tmp))
  {
    group=names(m.colors)[i]
    useCol=group
    if (group=="lightcyan") useCol="cadetblue1"
    if (group=="white") useCol="cadetblue1"
    if (group=="lightyellow") useCol="cadetblue1"
    boxplot(nm[dynamicColors==group,],outline=FALSE, border=useCol, names=names(nm),ylab="log2(Norm. TMT intensity)", main=paste("Cluster",i,'(','n=',length(tb[dynamicColors==group,1]),')',sep=''))
  }
  dev.off()
  
  # Calculate eigengenes
  MEList = moduleEigengenes(datExpr, colors = dynamicColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs);
  # Cluster module eigengenes
  
  #### I MODIFIED hclust TO CHECK=FALSE
  METree = hclust(as.dist(MEDiss), method = "average");
  
  # Call an automatic merging function
  MEDissThres = min_cluster_dist
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors;
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs;
  
  # Plot the result
  #pdf(file=paste(outputDir,'/eigenGene_tree.pdf',sep=''),width=10,height=10)
  #plot(METree, main = "Clustering of module eigengenes",
  #xlab = "", sub = "",check=FALSE)
  # Plot the cut line into the dendrogram
  #abline(h=MEDissThres, col = "red")
  #dev.off()
  
  
  
  # Check plothclust invalid dendrogram input #
  # Plot the dendrogram and colors underneath
  pdf(file=paste(outputDir,'/dendrogram_merged.pdf',sep=''),width=12,height=12)
  plotDendroAndColors(geneTree, mergedColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  dev.off()
  
  # hub genes
  # Calculate eigengenes again after merging
  MEList = moduleEigengenes(datExpr, colors = mergedColors)
  MEs = MEList$eigengenes
  #write.table(MEs,file='phosProteome_Cluster_trend_consensus_eigengenes.txt',quote=F,sep="\t")
  
  # kME
  datKME=signedKME(datExpr, MEs, outputColumnName="MM.")
  
  # attach kME to tb
  n=ncol(tb)
  tb[,(n+1):(n+ncol(datKME))]=datKME[,1:ncol(datKME)]
  
  #gn=whl$GN
  gn=whl[,2]
  tb$group=mergedColors
  tb$GN=gn
  
  # re-assign module membership by kME (>0.7)
  # question: how many genes with kME > 0.7 for more than one module?
  kME_cut=min_kME
  tboffset=ncol(whl)-col_offset
  k=rep(0,nrow(tb));
  for(i in 1:nrow(tb)) {
    for (j in (tboffset+1):(ncol(tb)-2)) {
      if (is.na(tb[i,j])==FALSE){
        if (tb[i,j]>=kME_cut) k[i]=k[i]+1
      }
    }
  }
  table(k)
  
  # strategy: try to find the 'best' module by assigning each gene to highest kME module with cutoff
  kME_cut=min_kME
  for(i in 1:nrow(tb)) {
    maxkme=max(tb[i,(tboffset+1):(ncol(tb)-2)])
    if (is.na(maxkme)==FALSE){
      if (maxkme>=kME_cut) {
        tb[i,ncol(tb)-1]=strsplit(colnames(tb)[which(tb[i,]==maxkme)],"\\.")[[1]][2]
      } else {
        tb[i,ncol(tb)-1]='grey'
      }
    }
  }
  
  # update mergedColors
  mergedColors.old=mergedColors
  mergedColors=tb[,ncol(tb)-1]
  
  # how many are changed by this additional step?
  k=0
  for (i in 1:length(mergedColors)) {
    if (mergedColors[i] != mergedColors.old[i]) k=k+1
  }
  k
  table(mergedColors)
  table(mergedColors.old)
  
  # Plot the dendrogram and colors underneath
  pdf(file=paste(outputDir,'/dendrogram_final.pdf',sep=''),width=12,height=12)
  plotDendroAndColors(geneTree, mergedColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  dev.off()
  
  # print genes (sorted by kME) for each module
  tmp=table(mergedColors);length(tmp)
  tmp=tmp[order(tmp,decreasing=T)]
  
  hasGrey=which(names(tmp)=='grey')
  if (length(hasGrey)>1) {
    group='grey'
    b=whl[mergedColors==group,1:ncol(whl)]
    write.table(b[!is.na(b[,1]),],file="notClustered.txt",quote=F,col.name=T,row.name=F,sep="\t")
  }
  
  if (length(hasGrey)>=1) tmp=tmp[-which(names(tmp)=='grey')] # remove 'grey' module
  
  m.colors=tmp
  tb1=whl
  
  for (i in 1:length(tmp))
  {
    group=names(m.colors)[i]
    colnm=paste("MM",group,sep='.');
    a=datKME[mergedColors==group,colnm]
    b=cbind(tb1[mergedColors==group,1:ncol(tb1)],a)
    write.table(b[order(a,decreasing = T),],file=paste(outputDir,'/',"Cluster_",i,".txt",sep=""),quote=F,col.name=T,row.name=F,sep="\t")
  }
  
  # print one file of WPC information
  tb1$cluster='n/a'
  for (i in 1:length(tmp)) {
    group=names(m.colors)[i]
    tb1$cluster[mergedColors==group]=i
  }
  filename1 <- sprintf('%s/Coexpression_clustered_norm_exp_matrix.txt', outputDir)
  write.table(tb1[,c(ncol(tb1),1:(ncol(tb1)-1))],file=filename1,quote=F,col.name=T,row.name=F,sep="\t")
  
  # print out correlation table (vs. kME) for each gene
  filename2 <- sprintf('%s/Coexpression_clustered_norm_exp_matrix_kMEcor.txt', outputDir)
  write.table(cbind(tb1[,c(ncol(tb1),1:(ncol(tb1)-1))],tb[,grep('MM',colnames(tb),fixed=T)]),file=filename2,quote=F,col.name=T,row.name=F,sep="\t")
  
  # print boxplot
  
  # normalize signal to print trend (like what STEM do)
  nm=tb[,1:plex]
  for (i in 1:plex) {
    nm[,i]=tb[,i]-tb[,1]
  }
  summary(nm)
  head(nm)
  dim(nm)
  
  # added on 1/23/19
  # if too many (>25) modules generated, only print the top 25
  if (length(tmp)>25) tmp=tmp[1:25]
  # added end
  
  png(file=paste(outputDir,'/boxplot_coexpression_cluster.png',sep=''),width=500,height=500)
  box_col=ceiling(sqrt(length(tmp)))
  box_row=ceiling(length(tmp)/box_col)
  par(las=2,mfrow=c(box_row,box_col),mar = c(1, 1, 1, 1),cex.lab=1.2)
  for (i in 1:length(tmp))
  {
    group=names(m.colors)[i]
    useCol=group
    if (group=="lightcyan") useCol="cadetblue1"
    if (group=="white") useCol="cadetblue1"
    if (group=="lightyellow") useCol="cadetblue1"
    boxplot(nm[mergedColors==group,],outline=FALSE,border =useCol, names=names(nm),ylab="log2(Norm. TMT intensity)", main=paste("Cluster",i,'(','n=',length(tb[mergedColors==group,1]),')',sep=''))
  }
  dev.off()
  
  for (i in 1:length(tmp)) {
    group=names(m.colors)[i]
    useCol=group
    if (group=="lightcyan") useCol="cadetblue1"
    if (group=="white") useCol="cadetblue1"
    if (group=="lightyellow") useCol="cadetblue1"
    
    # small size
    png(file=paste(outputDir,'/small_boxplot_coexpression_cluster_',i,'.png',sep=''),width=500,height=500)
    par(las=2,mfrow=c(1,1),mar = c(1, 1, 1, 1))
    
    boxplot(nm[mergedColors==group,],outline=FALSE,border =useCol, names=names(nm),ylab="log2(Norm. TMT intensity)", main=paste("Cluster",i,'(','n=',length(tb[mergedColors==group,1]),')',sep=''))
    dev.off()
    
    # large size
    png(file=paste(outputDir,'/large_boxplot_coexpression_cluster_',i,'_',group,'.png',sep=''),width=1000,height=1000)
    par(las=2,mfrow=c(1,1),mar = c(1, 1, 1, 1))
    
    boxplot(nm[mergedColors==group,],outline=FALSE,border =useCol, names=names(nm),ylab="log2(Norm. TMT intensity)", main=paste("Cluster",i,'(','n=',length(tb[mergedColors==group,1]),')',sep=''))
    dev.off()
  }
  
  # print pretty trends for each cluster
  # verage replicates first
  # for each protein
  print(nrow(datExpr))
  datExpr2=matrix(nrow=nrow(datExpr)/2,ncol=ncol(datExpr))
  colnames(datExpr2)=colnames(datExpr)
  for (i in 1:ncol(datExpr2)) {
    for (j in 1:nrow(datExpr2)) {
      datExpr2[j,i]=mean(datExpr[(j*2-1):(j*2),i])
    }
  }
  if (nrow(datExpr)==4 || nrow(datExpr)==5){
    datExpr2=matrix(nrow=nrow(datExpr),ncol=ncol(datExpr))
    colnames(datExpr2)=colnames(datExpr)
    for (i in 1:ncol(datExpr2)) {
      for (j in 1:nrow(datExpr2)) {
        datExpr2[j,i]=datExpr[j,i]
      }
    }
  }
  dim(datExpr2)
  
  
  # eigenGene
  # Calculate eigengenes again after re-assign
  MEList = moduleEigengenes(datExpr, colors = mergedColors)
  MEs = MEList$eigengenes
  MEs2=matrix(nrow=nrow(datExpr2),ncol=ncol(MEs))
  rownames(MEs2)=rownames(datExpr2)
  colnames(MEs2)=colnames(MEs)
  for (i in 1:ncol(MEs2)) {
    for (j in 1:nrow(datExpr2)) {
      MEs2[j,i]=mean(MEs[(j*2-1):(j*2),i])
    }
  }
  MEs2
  # step 1: Data normalization (mean = 0, SD = 1):
  # for each protein 
  tbn=scale(datExpr2)
  tbn=t(tbn)
  # For each eigenGene
  MEn=scale(MEs2)
  # step 2: Configure color:
  # range: 0.7-1
  # color: yellow . green . blue - red
  library(gplots)
  my_palette <- colorRampPalette(c("green", "blue","red"))(n = 301)
  # step3: Plot trend:
  length(tmp)
  #skipRange=c(2,3)# for plot that combines WT and KO, like Fig. 5B in Tan et al., Immunity 2017
  png(file=paste(outputDir,'/trends_coexpression_cluster.png',sep=''),width=500,height=500)
  box_col=ceiling(sqrt(length(tmp)))
  box_row=ceiling(length(tmp)/box_col)
  par(las=2,mfrow=c(box_row,box_col),mar = c(2, 2, 2, 2),cex.lab=.8)
  # for each cluster:
  for (i in 1:length(tmp))
  {
    print('started')
    print(length(tmp))
    group=names(m.colors)[i]
    print(group)
    m=tbn[mergedColors==group,]
    
    # get kME
    colnm=paste("MM",group,sep='.');
    colcn=grep(colnm,colnames(datKME),fixed=T)
    a=datKME[mergedColors==group,colcn]
    print(a)
    print(typeof(a))
    if (is.double(a)==FALSE){
      a <- as.double(a[,1])
    }
    # plot the framework
    colnm=paste("ME",group,sep='');
    colcn=grep(colnm,colnames(MEn),fixed=T)
    if (length(colcn)>1){
      colcn <- colcn[1]
    }
    specific_MEn <- MEn[!is.na(MEn[,colcn]),colcn]
    plot(MEn[,colcn],type='l',ylab='Relative level',xlab='',ylim=c(min(specific_MEn)-0.5,max(specific_MEn)+0.5),main=paste("Cluster",i,' (n=',length(tb[mergedColors==group,1]),')',sep=''))
    #lines(skipRange,MEn[skipRange,colcn],col='white')
    
    # for each protein, 
    print(m)
    print(typeof(m))
    for (j in 1:nrow(m)) {
      # color determined by membership
      crtCol=my_palette[round((a[j] - 0.7)*1000)]
      # plot trend
      lines(m[j,],col=crtCol)
      #lines(skipRange,m[j,skipRange],col='white')
    }
    # plot consensus trend (eigenGene)
    # color: black
    #lines(MEn[,colcn],type='l')
  }
  dev.off()
  
  for (i in 1:length(tmp))
  {
    png(file=paste(outputDir,'/trend_coexpression_cluster_',i,'.png',sep=''),width=500,height=500)
    par(mar = c(5, 2, .75, 1),cex.lab=1.2)
    group=names(m.colors)[i]
    m=tbn[mergedColors==group,]

    # get kME
    colnm=paste("MM",group,sep='.');
    colcn=grep(colnm,colnames(datKME),fixed=T)
    a=datKME[mergedColors==group,colcn]
    if (is.double(a)==FALSE){
      a <- as.double(a[,1])
    }
    # plot the framework
    colnm=paste("ME",group,sep='');
    colcn=grep(colnm,colnames(MEn),fixed=T)
    if (length(colcn)>1){
      colcn <- colcn[1]
    }
    specific_MEn <- MEn[!is.na(MEn[,colcn]),colcn]
    plot(MEn[,colcn],type='l',ylab='Relative level',xlab='',ylim=c(min(specific_MEn)-0.5,max(specific_MEn)+0.5),main=paste("Cluster",i,' (n=',length(tb[mergedColors==group,1]),')',sep=''))
    #lines(skipRange,MEn[skipRange,colcn],col='white')
    
    # for each protein, 
    for (j in 1:nrow(m)) {
      # color determined by membership
      crtCol=my_palette[round((a[j] - 0.7)*1000)]
      # plot trend
      lines(m[j,],col=crtCol)
      #lines(skipRange,m[j,skipRange],col='white')
    }
    # plot consensus trend (eigenGene)
    # color: black
    #lines(MEn[,colcn],type='l')
    dev.off()
  }
  
  # make a new matrix for eigen genes
  mes3 = matrix(nrow=nrow(MEs),ncol=length(tmp))
  colnames(mes3)=names(tmp)
  rownames(mes3)=rownames(datExpr)
  
  for (i in 1:length(tmp)) {
    group=names(m.colors)[i]
    colnm=paste("ME",group,sep='');
    colcn=grep(colnm,colnames(MEn),fixed=T)
    mes3[,i]=MEs[,colcn]
  }
  
  write.table(mes3,file=paste(outputDir,"/MEs_final.txt", sep=''),quote=F,col.name=T,row.name=T,sep="\t")
  
  # print barplot for eigen genes
  pdf(paste(outputDir,'/barplot_eigengene.pdf',sep=''),15,15)
  for (k in 1:ncol(mes3)) {
    par(las=3,mar=c(16,4,4,1))
    barplot(mes3[,k],main=paste('cluster ',k,' (n=',tmp[k],')',sep=''),col=names(tmp)[k])
  }
  dev.off()
  
  #--------------------
  # parse meta data
  print(metaFile)
  if (is.null(metaFile)==FALSE) {
    print('Check1Passed')
    metaFile <- metaFile$datapath
    extension <- file_ext(metaFile)
    if (extension == 'xlsx'){                                                                            ### If the User File is .xlsx ...
      grp <- read_excel(metaFile, na='NA')                                                               ### Call read_excel function on the user input file
    } else if (extension == 'csv'){                                                                     ### If the User file is .csv  ...
      grp <- read.csv(metaFile)                                                               ### Also call read_excel function
    } else if (extension == 'txt'){                                                                     ### If the User file is ,txt
      grp=read.table(metaFile,head=T,sep="\t") # rownames should match sample names      
    }        
    summary(grp)
    head(grp)
    dim(grp)
    
    #----- prepare assayData and pData
    # meta data
    pData = grp
    rownames(pData)=grp[,1]
    # assay data
    assayData=t(mes3)
    print(sprintf('nrow mes3 %d nrow pdata %d', nrow(mes3),nrow(pData)))
    print(rownames(mes3))
    print(rownames(pData))
    new_rownames <- c()
    rowname_vector <- as.vector(rownames(mes3))
    for (sample in 1:length(rowname_vector)){
      new_rownames <- c(new_rownames,substr(rowname_vector[sample],3,nchar(rowname_vector[sample])-1))
    }
    print(new_rownames)
    print(sum(new_rownames==rownames(pData)))
    # check sample name consistency
    if (sum(new_rownames==rownames(pData))==nrow(pData)) {
      print('Check2Passed')
      #-------------------- boxpot: each cluster vs. each factor
      # 1 file 1 factor; for each factor: 1 cluster per page
      egM=assayData
      colUsed=2:ncol(pData) # use all columns (excepct the 1st one, which is the sample IS)
      print('Length of Cols Used')
      print(length(colUsed))
      for (i in 1:length(colUsed)) { # for each factor: 1 file 1 factor
        #if (is.factor(pData[,colUsed[i]])) {
          print('Check3Passed')
          png(paste(outputDir,'/eigengene_vs_',colnames(pData)[colUsed[i]],'.png',sep=''),1000,800)
          box_e_col=ceiling(sqrt(nrow(egM)))
          box_e_row=ceiling(nrow(egM)/box_col)
          par(las=2,mfrow=c(box_e_row,box_e_col),mar=c(6,6,6,6),cex.main=3,cex.lab=2,cex.axis=2)
          for (k in 1:nrow(egM)) { # for each factor: 1 cluster per page
            #par(las=3,mar=c(1,1,1,1))
            par(cex=2)
            boxplot(egM[k,]~pData[,colUsed[i]],main=paste('cluster ',k,' (n=',tmp[k],')',sep=''), xlab=colnames(pData)[colUsed[i]], ylab='Relative Abundance')
            stripchart(egM[k,]~pData[,colUsed[i]], vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue')
          }
          dev.off()
        #}
      }
    } else {
        warning("meta data sample IDs does NOT match the provided quantification matrix!\n(eigen analysis will be skipped)")
    }
  }
  #save.image(file=paste(outputDir,'/wgcna2.RData', sep=''))
}
