#################
### Libraries ###
#################

library("Seurat")
library("SingleCellExperiment")
library("scran")
library("scater")
library("bluster")

##################
### Data Input ###
##################

# read files
filename = file.choose()
hd = readRDS(filename)
filename = file.choose()
active = readRDS(filename)
filename = file.choose()
recovered = readRDS(filename)

# Feature Selection

active <- FindVariableFeatures(active, selection.method = "vst", nfeatures = 2000)
recovered <- FindVariableFeatures(recovered, selection.method = "vst", nfeatures = 2000)
hd <- FindVariableFeatures(hd, selection.method = "vst", nfeatures = 2000)

top10=head(VariableFeatures(active),10)
plot1=VariableFeaturePlot(active)
LabelPoints(plot=plot1,points=top10,repel=TRUE)
#scaling

all.genes_active <- rownames(active)
active <- ScaleData(active, features = all.genes_active)

all.genes_recovered <- rownames(recovered)
recovered <- ScaleData(recovered, features = all.genes_recovered)

all.genes_hd <- rownames(hd)
hd <- ScaleData(hd, features = all.genes_hd)

#PCA

active <- RunPCA(active, features = VariableFeatures(object = active))

recovered <- RunPCA(recovered, features = VariableFeatures(object = recovered))

hd <- RunPCA(hd, features = VariableFeatures(object = hd))

#Clustering

active <- FindNeighbors(active, dims = 1:10)
active <- FindClusters(active, resolution = 0.5)

recovered <- FindNeighbors(recovered, dims = 1:10)
recovered <- FindClusters(recovered, resolution = 0.5)

hd <- FindNeighbors(hd, dims = 1:10)
hd <- FindClusters(hd, resolution = 0.5)

#UMAP

active <- RunUMAP(active, dims = 1:10)

recovered <- RunUMAP(recovered, dims = 1:10)

hd <- RunUMAP(hd, dims = 1:10)

#Visualisierung von UMAp

DimPlot(active, reduction = "umap")

DimPlot(recovered, reduction = "umap")

DimPlot(hd, reduction = "umap")

####################################################################################
#####Finding differentially expressed features (cluster biomarkers)#################
####################################################################################

#FindMarkers in PB Cluster

##active##

Cluster0.active <- FindMarkers(active, ident.1 = 0, min.pct = 0.25)
Cluster1.active <- FindMarkers(active, ident.1 = 1, min.pct = 0.25)
Cluster2.active <- FindMarkers(active, ident.1 = 2, min.pct = 0.25)
Cluster3.active <- FindMarkers(active, ident.1 = 3, min.pct = 0.25)
Cluster4.active <- FindMarkers(active, ident.1 = 4, min.pct = 0.25)
Cluster5.active <- FindMarkers(active, ident.1 = 5, min.pct = 0.25)
Cluster6.active <- FindMarkers(active, ident.1 = 6, min.pct = 0.25)

##recovered##

Cluster0.recovered <- FindMarkers(recovered, ident.1 = 0, min.pct = 0.25)
Cluster1.recovered <- FindMarkers(recovered, ident.1 = 1, min.pct = 0.25)
Cluster2.recovered <- FindMarkers(recovered, ident.1 = 2, min.pct = 0.25)
Cluster3.recovered <- FindMarkers(recovered, ident.1 = 3, min.pct = 0.25)
Cluster4.recovered <- FindMarkers(recovered, ident.1 = 4, min.pct = 0.25)
Cluster5.recovered <- FindMarkers(recovered, ident.1 = 5, min.pct = 0.25)
Cluster6.recovered <- FindMarkers(recovered, ident.1 = 6, min.pct = 0.25)

##hd##

Cluster0.hd <- FindMarkers(hd, ident.1 = 0, min.pct = 0.25)
Cluster1.hd <- FindMarkers(hd, ident.1 = 1, min.pct = 0.25)
Cluster2.hd <- FindMarkers(hd, ident.1 = 2, min.pct = 0.25)

#FeaturePlot

FeaturePlot(active, features = c("IGHD","IGHM","IGHA2","IGHA1","IGHG1","IGHG2","BACH2","EBF-1"))


#####################################################################################################
#####Bis hier funktioniert allles super. SingleR macht dann Probleme###############################
# Als Referenz habe ich 4 Datensets ausprobiert
#StoeckiusHashingData(mode='human') aus dem scRNAseq package
#KotliarovPBMCData() aus dem scRNAseq Package
#MairPBMCData() aus dem scRNAseq Package
#HumanPrimaryCellAtlasData() aus celldex package
#DatabaseImmuneCellExpressionData() aus celldex
# keins der Datensets funktionierte und jedes Set brachte eigene Probleme und Fehlermeldungen.

#SingleR Celltype Annotation##

#library laden
library(SingleR)
library(celldex)
library(scRNAseq)


#Referenz definieren

ref=celldex::HumanPrimaryCellAtlasData()
ref=logNormCounts(ref)

#Results

Label_active= SingleR(test=active,ref=ref,labels=ref$label.main)

