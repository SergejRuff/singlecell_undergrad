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
  
  
  #########################################
  ### Differentially Expressed Features ###
  #########################################
  
  clusterList = vector(mode = "list", length = length(levels(SeuratObj$seurat_clusters)))
  
  #Loop that finds markers, then stores them in an additional cluster list 
  for(i in 1:length(clusterList)) clusterList[[i]] = FindMarkers(SeuratObj, ident.1 = as.numeric(levels(SeuratObj$seurat_clusters)[i]), min.pct = 0.25)
  
  
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


##################
### Data Input ###
##################

# read files
active = readRDS("A:/RFILES/E-MTAB-11011 - Single Cell Analysis of B cells in COVID-19 comprising active and recovered disease/.rds-files der Patienten/active-Covid-Kranke/pbmc.active.2.5.3.8_gex_and_vdj.rds")


########################################
### Remove wrongfully assigned cells ###
########################################

hpca.se = celldex::HumanPrimaryCellAtlasData() #reference data
active = remWC(active)


########################
### Feature Analysis ###
########################


active_data = scsPipeline(active)
active = active_data[[1]]
active_clusters = active_data[[2]]

#################################################
#########Daten als csv exportieren#####
#################################################

setwd("A:/RFILES/FindMarkers")

cluster_list=as.vector(active_clusters,mode="list")

for (i in 1:length(cluster_list)){
  cluster_list[[i]]=cbind(cluster_list[[i]],i)
  cluster_list[[i]]=cbind(cluster_list[[i]],rownames(cluster_list[[i]]))
}

cluster_list=do.call(rbind,cluster_list)

write.table(cluster_list,"DESeq2ErrorFindMarkers.csv",quote=FALSE,sep=";")

##############################################################################################################
##############################################################
################csv-analyse###################################
##############################################################

boot=cluster_list

##############Anzahl an Genen und unique Genen vergleichen#####################
nrow(boot)

length(unique(boot$`rownames(cluster_list[[i]])`))
