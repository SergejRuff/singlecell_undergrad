### This code follows a basic pipeline for Single Cell Sequencing Data Analysis
### based on the example dataset E-MTAB-11011

setwd("A:/RFILES/setwd")


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
  
 
  
  ##################
  ### Clustering ###
  ##################
  
  SeuratObj = FindNeighbors(SeuratObj, dims = 1:10)
  SeuratObj = FindClusters(SeuratObj, resolution = 0.5)
  

  #########################################
  ### Differentially Expressed Features ###
  #########################################
  
  clusterList = vector(mode = "list", length = length(levels(SeuratObj$seurat_clusters)))
  
  #Loop that finds markers, then stores them in an additional cluster list 
  for(i in 1:length(clusterList)) clusterList[[i]] = FindMarkers(SeuratObj, ident.1 = as.numeric(levels(SeuratObj$seurat_clusters)[i]), min.pct = 0.25)
   

 
  
  #To return multiple objects, wrap them in a list first
  return(clusterList)
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

########################################
###########Bootstrap####################
########################################
#as.matrix(GetAssayData(active))
#n=round(ncol(A)/100*80)
len=100
iteration=vector(mode="list",length=len)
n=round(ncol(active)/100*80)
start=Sys.time()
for(i in 1:len){
  s=sample(1:ncol(active),size=n)
  B=active[,s]
  iteration[[i]]=scsPipeline(B)
}
end=Sys.time()
print(end-start)
#################################################
#########Bootstrap-Daten als csv exportieren#####
#################################################

setwd("A:/RFILES/bootstrap/15062022")

for (i in 1:length(iteration)){
  for(j in 1:length(iteration[[i]])){
    iteration[[i]][[j]]=cbind(iteration[[i]][[j]],i,j)
    iteration[[i]][[j]]=cbind(iteration[[i]][[j]],rownames(iteration[[i]][[j]]))
  }
}


for(i in 1:length(iteration))iteration[[i]]= do.call(rbind,iteration[[i]])
result=do.call(rbind,iteration)
write.table(result,"100Bootstraping150622.csv",quote=FALSE,sep=";")


  