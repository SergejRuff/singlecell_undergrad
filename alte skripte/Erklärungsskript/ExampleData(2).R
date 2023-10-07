
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

# merge all data
#merged = merge(hd, y = c(active, recovered), add.cell.ids = c("hd", "active", "recovered"))
#merged = merge(active, y = recovered,add.cell.ids = c("active", "recovered"))

# create Single Cell Experiment class
#sce = as.SingleCellExperiment(merged)

sce_active = as.SingleCellExperiment(active)
sce_recovered = as.SingleCellExperiment(recovered)

#sce = sce_active ...

#######################
### Quality Control ###
#######################

# identify cells with more than 10% mitochondrial DNA content and remove those cells
mt = which(sce$percent.mt >= 10)
if(length(mt) > 0) sce = sce[,-mt] 

# identify cells with more than 2500 genes or less than 200 genes
gpc = colSums(logcounts(sce) > 0)
gpcr = which(gpc > 2500 | gpc < 200)
if(length(gpcr) > 0) sce = sce[,-gpcr]

###########################
### Data Pre-Processing ###
###########################

# log-normalize counts
#sce = logNormCounts(sce)

################################################################################

################
### Analyses ###
################

# Feature Selection
dec = modelGeneVar(sce)
hvg = getTopHVGs(dec, prop=0.1)

# Dimension Reduction
sce = runPCA(sce, ncomponents = 10)

# Clustering
set.seed(666)
colLabels(sce) = clusterCells(sce, use.dimred='PCA', BLUSPARAM=NNGraphParam(cluster.fun="louvain")) 

# Visualisation
sce = runUMAP(sce, dimred = 'PCA')
plotUMAP(sce, colour_by="label")

# Find Markers
markers = findMarkers(sce, test.type="wilcox", direction="up", lfc=1)