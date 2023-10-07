#SingleCell RNA-seq practice

##############
#packages#####
##############

library("Seurat")


###########################################################
############Seurat Analysis################################
###########################################################


##############
#Import File##
##############

# this shows how to import mtx data: you need features and barcodes as seperate files and a matrix-file (mtx)
Sheep1<-ReadMtx(mtx="A:/RFILES/Nematoda/schaf1/matrix.mtx",features="A:/RFILES/Nematoda/schaf1/features/features.tsv",cells="A:/RFILES/Nematoda/schaf1/barcode/barcodes.tsv")
Sheep1<-CreateSeuratObject(counts=Sheep1)

Sheep2<-ReadMtx(mtx="A:/RFILES/Nematoda/schaf2/matrix.mtx",features="A:/RFILES/Nematoda/schaf2/features/features.tsv",cells="A:/RFILES/Nematoda/schaf2/barcode/barcodes.tsv")
Sheep2<-CreateSeuratObject(counts=Sheep2)


##########################################
#Quality-Control##########################
##########################################

#calculate percantage of mitochndirial RNA and add new Col.
Sheep1[["percent.mt"]] <- PercentageFeatureSet(Sheep1, pattern = "^MT-")
Sheep2[["percent.mt"]] <- PercentageFeatureSet(Sheep2, pattern = "^MT-")

Sheep1<-subset(Sheep1, subset = nFeature_RNA > 300 & nFeature_RNA < 6500 & percent.mt < 20)
Sheep2<-subset(Sheep2, subset = nFeature_RNA > 300 & nFeature_RNA < 5750 & percent.mt < 20)

######################################
######Normalize#######################
######################################

Sheep1 <- NormalizeData(Sheep1, normalization.method = "LogNormalize", scale.factor = 10000)
Sheep2 <- NormalizeData(Sheep2, normalization.method = "LogNormalize", scale.factor = 10000)

Sheep1 <- NormalizeData(Sheep1)
Sheep2 <- NormalizeData(Sheep2)