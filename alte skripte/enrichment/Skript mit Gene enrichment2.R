### This code follows a basic pipeline for Single Cell Sequencing Data Analysis
### based on the example dataset E-MTAB-11011

setwd("A:/RFILES/für 25.05.22")


#################
### Libraries ###
#################

library("Seurat")
library("SingleCellExperiment")
library("scran")
library("scater")
library("bluster")
library("celldex")
library("SingleR")
library("dplyr")
library("biomaRt")
library("clusterProfiler")
library("org.Hs.eg.db")

#################
### Functions ###
#################

# function to compute all steps in a basic single cell experiment pipeline
scsPipeline = function(SeuratObj){
  
  #########################
  ### Feature Selection ###
  #########################
  
  SeuratObj = FindVariableFeatures(SeuratObj, selection.method = "vst", nfeatures = 2000)

  ####################
  ### Data Scaling ###
  ####################
  
  all.genes_SeuratObj = rownames(SeuratObj)
  SeuratObj = ScaleData(SeuratObj, features = all.genes_SeuratObj)
  
  #################################
  ### Dimension Reduction (PCA) ###
  #################################
  
  SeuratObj = RunPCA(SeuratObj, features = VariableFeatures(object = SeuratObj))
  
  ##################
  ### Clustering ###
  ##################
  
  SeuratObj = FindNeighbors(SeuratObj, dims = 1:10)
  SeuratObj = FindClusters(SeuratObj, resolution = 0.5)
  
  ##################################
  ### Dimension Reduction (UMAP) ###
  ##################################
  
  SeuratObj = RunUMAP(SeuratObj, dims = 1:10)
  
  #####################
  ### Visualisation ###
  #####################
  
  dp=DimPlot(SeuratObj, reduction = "umap")
  print(dp)
  #########################################
  ### Differentially Expressed Features ###
  #########################################
  
  clusterList = vector(mode = "list", length = length(levels(SeuratObj$seurat_clusters)))
  
  #Loop that finds markers, then stores them in an additional cluster list 
  for(i in 1:length(clusterList)) clusterList[[i]] = FindMarkers(SeuratObj, ident.1 = as.numeric(levels(SeuratObj$seurat_clusters)[i]), min.pct = 0.25)
   

  #FeaturePlot
  fp=FeaturePlot(SeuratObj, features = c("IGHD","IGHM","IGHA2","IGHA1","IGHG1","IGHG2","BACH2","EBF1"), cols = c("cadetblue3", "black"))
  print(fp)
  
  #To return multiple objects, wrap them in a list first
  return(list(SeuratObj, clusterList))
}

#function to remove wrongly assigned cells from a single cell experiment stored as Seurat Object
remWC = function(SeuratObj, regex = ".*B_cell.*"){
  
  ############################
  ### Cell Type Annotation ###
  ############################
  
  #active 
  sce = as.SingleCellExperiment(SeuratObj) #convert Seurat to SingleCellExperiment
  se = as(sce, "SummarizedExperiment") #convert SingleCellExperiment  to SummarizedExperiment
  rownames(se) = rownames(sce) #re-add the rownames 
  
  # predict cell type
  pred = SingleR(test = se, ref = hpca.se, assay.type.test=1,
                 labels = hpca.se$label.main)
  
  pred$first.labels #show results
  
  #percentage b cells
  pB = sum(grepl(".*B_cell.*", pred$first.labels)) / length(pred$first.labels) * 100
  print(paste0("Percentage B-cells: ", round(pB, 2), " %"))
  
  SeuratObj$CellAnnotation = pred$first.labels
  
  ###################################
  ### Removal of wrong cell types ###
  ###################################
  
  ca = as.vector(SeuratObj$CellAnnotation) #vector with cell annotations
  indices = which(grepl(regex, ca)) #indices of cells to keep
  SeuratObj = SeuratObj[,indices] #retain cells according to indices
  
  return(SeuratObj) #return Seurat Object after removal of wrongly assigned cells
}

#function to export the top200 genes for each cluster 
exT200 = function(clusterList, name, path = getwd()){
  for(i in 1:length(clusterList)){
    #take the row names from the first 200 rows in clusters unless there's fewer rows available
    num = ifelse(nrow(clusterList[[i]]) < 200, nrow(clusterList[[i]]), 200)
    top200 = rownames(clusterList[[i]])[1:num]
    write.table(top200, file = paste0(path,"/", name, "_", i, "_top200.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}

################################################
########Gene enrichment(KEGG)###################
################################################

#Gene Set Enrichment Analysis based on log-fold change expression values using KEGG
GE = function(GeneList, pval = 1, byNames = FALSE, n = 200){

  GeneList = cbind(rownames(GeneList), GeneList) #append the row names (genes) as column
  colnames(GeneList)[1] = "Genename"
  
  # find corresponding entrez ID for genes
  genes = mapIds(org.Hs.eg.db, keys = GeneList[,1],
                 column = "ENTREZID",
                 keytype = "SYMBOL") 
  
  # create data frame with gene name and entrez id as columns
  genes = as.data.frame(cbind(names(genes), as.vector(genes))) 
  colnames(genes) = c("Genename", "entrezID")
  Merged_data = left_join(GeneList, genes, by = "Genename") #combine tables
  
  if(byNames == TRUE){ #find gene enrichments by names only
    Merged_data = Merged_data[order(Merged_data$p_val_adj, decreasing = FALSE),] #sort by increasing adj. p value
    num = ifelse(nrow(Merged_data) < n, nrow(Merged_data), n) #take the first n genes or all genes if num < n
    topGenes = Merged_data$entrezID[1:num] #retain the first n genes 
    
    kegg = enrichKEGG(topGenes,
                      organism = "hsa",
                      pvalueCutoff = pval)
  } else { 
    #find gene enrichments based on all genes and their log-fold change expression
    ngl = Merged_data$avg_log2FC #extract the log-fold change in expression
    
    names(ngl) = Merged_data$entrezID
    
    ngl = ngl[-which(is.na(names(ngl)))] #remove genes without entrez ID
    ngl = sort(ngl, decreasing = TRUE) #sort by increasing fold change value
    
    #perform enrichment analysis
    kegg = gseKEGG(geneList= ngl,
                   organism= "hsa",
                   minGSSize= 5,
                   pvalueCutoff= pval,
                   verbose= TRUE)
  }
  
  return(data.frame(kegg))
}


##################
### Data Input ###
##################

# read files
hd = readRDS("A:/RFILES/E-MTAB-11011 - Single Cell Analysis of B cells in COVID-19 comprising active and recovered disease/.rds-files der Patienten/HD-GesundePatienten/pbmc.HD_gex_and_vdj.rds")
active = readRDS("A:/RFILES/E-MTAB-11011 - Single Cell Analysis of B cells in COVID-19 comprising active and recovered disease/.rds-files der Patienten/active-Covid-Kranke/pbmc.active.2.5.3.8_gex_and_vdj.rds")
recovered = readRDS("A:/RFILES/E-MTAB-11011 - Single Cell Analysis of B cells in COVID-19 comprising active and recovered disease/.rds-files der Patienten/recovered-Erholte Patienten/pbmc.recovered.14.16.26_gex_and_vdj.rds")

########################################
### Remove wrongfully assigned cells ###
########################################

hpca.se = celldex::HumanPrimaryCellAtlasData() #reference data
active = remWC(active)
recovered = remWC(recovered)
hd = remWC(hd)

########################
### Feature Analysis ###
########################

hd_data = scsPipeline(hd)
hd = hd_data[[1]]
hd_clusters = hd_data[[2]]

active_data = scsPipeline(active)
active = active_data[[1]]
active_clusters = active_data[[2]]

recovered_data = scsPipeline(recovered)
recovered = recovered_data[[1]]
recovered_clusters = recovered_data[[2]]

####################################
### Gene Set Enrichment Analysis ###
####################################

# GSEA for each cluster based on log-fold change in expression using KEGG

#for each list item (cluster) perform a separate gene set enrichment analysis
active_kegg = lapply(active_clusters, function(x) GE(x))
recovered_kegg = lapply(recovered_clusters, function(x) GE(x))
hd_kegg = lapply(hd_clusters, function(x) GE(x))

names(active_kegg) = names(active_clusters)
names(recovered_kegg) = names(recovered_clusters)
names(hd_kegg) = names(hd_clusters)


######################################
### Gene Set Enrichment Analysis 2 ###
######################################

# GSEA for each cluster based on top 200 differentially expressed genes using KEGG 

active_kegg2 = lapply(active_clusters, function(x) GE(x, byNames = TRUE, n = 200))
recovered_kegg2 = lapply(recovered_clusters, function(x) GE(x, byNames = TRUE, n = 200))
hd_kegg2 = lapply(hd_clusters, function(x) GE(x, byNames = TRUE, n = 200))

names(active_kegg2) = names(active_clusters)
names(recovered_kegg2) = names(recovered_clusters)
names(hd_kegg2) = names(hd_clusters)




top 200 genes for each cluster
exT200(hd_clusters, "hd")
exT200(active_clusters, "active")
exT200(recovered_clusters, "recovered")


  
  