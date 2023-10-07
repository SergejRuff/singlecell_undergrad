### This code follows a basic pipeline for Single Cell Sequencing Data Analysis

setwd("A:/RFILES/für Mittwoch 18.05.2022")

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
  
  #Loop that finds markers, then stores them in the Seurat Object as Metadata 
  for(i in 1:length(clusterList)) clusterList[[i]] = FindMarkers(SeuratObj, ident.1 = as.numeric(levels(SeuratObj$seurat_clusters)[i]), min.pct = 0.25)
   

  #FeaturePlot
  fp=FeaturePlot(SeuratObj, features = c("IGHD","IGHM","IGHA2","IGHA1","IGHG1","IGHG2","BACH2","EBF1"), cols = c("cadetblue3", "black"))
  print(fp)
  
  ############################
  ### Cell Type Annotation ###
  ############################
  
  hpca.se = celldex::HumanPrimaryCellAtlasData() #reference data
  
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
  
  SeuratObj$CellAnnotation = pred$first.labels
  
  #To return multiple objects, wrap them in a list first
  list = vector(mode = "list", length = 2) # list with two entries
  list[[1]] = SeuratObj # Seurat Object
  list[[2]] = clusterList # List with Clusters for the corresponding Seurat Object 
  
  return(list)
}

#function to remove wrongly assigned cells from a single cell experiment stored as Seurat Object

remWC = function(SeuratObj, regex = ".*B_cell.*"){
  
  ca = as.vector(SeuratObj$CellAnnotation) #vector with cell annotations
  indices = which(grepl(regex, ca)) #indices of cells to keep
  SeuratObj = SeuratObj[,indices] #retain cells according to indices
  
  return(SeuratObj) #return Seurat Object after removal of wrongly assigned cells
}

#function to export the top200 genes for each cluster 
exT200 = function(clusterList, name, path = getwd()){
  for(i in 1:length(clusterList)){
    #take the row names from the first 200 rows in hd_clusters unless there's fewer rows available
    num = ifelse(nrow(clusterList[[i]]) < 200, nrow(clusterList[[i]]), 200)
    top200 = rownames(clusterList[[i]])[1:num]
    write.table(top200, file = paste0(path,"/", name, "_", i, "_top200.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}

##################
### Data Input ###
##################

# read files
hd = readRDS("A:/RFILES/E-MTAB-11011 - Single Cell Analysis of B cells in COVID-19 comprising active and recovered disease/.rds-files der Patienten/HD-GesundePatienten/pbmc.HD_gex_and_vdj.rds")
active = readRDS("A:/RFILES/E-MTAB-11011 - Single Cell Analysis of B cells in COVID-19 comprising active and recovered disease/.rds-files der Patienten/active-Covid-Kranke/pbmc.active.2.5.3.8_gex_and_vdj.rds")
recovered = readRDS("A:/RFILES/E-MTAB-11011 - Single Cell Analysis of B cells in COVID-19 comprising active and recovered disease/.rds-files der Patienten/recovered-Erholte Patienten/pbmc.recovered.14.16.26_gex_and_vdj.rds")

hd_data = scsPipeline(hd)
hd = hd_data[[1]]
hd_clusters = hd_data[[2]]

active_data = scsPipeline(active)
active = active_data[[1]]
active_clusters = active_data[[2]]

recovered_data = scsPipeline(recovered)
recovered = recovered_data[[1]]
recovered_clusters = recovered_data[[2]]

########################################
### Remove wrongfully assigned cells ###
########################################

active = remWC(active)
recovered = remWC(recovered)
hd = remWC(hd)

###################
### Data Export ###
###################

#top 200 genes for each cluster
exT200(hd_clusters, "hd")
exT200(active_clusters, "active")
exT200(recovered_clusters, "recovered")