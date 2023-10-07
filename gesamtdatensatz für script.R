### This code follows a basic pipeline for Single Cell Sequencing Data Analysis
### based on the example dataset E-MTAB-11011

setwd("A:/RFILES/für 25.05.22")


#################
### Libraries ###
#################

library("Seurat")
library("celldex")
library("SingleR")
library("data.table")

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
########################################################################################################
#####################################
####filtering dataset ##########
#####################################

ClusterMarkers=as.vector(active_clusters,mode="list")

for (i in 1:length(ClusterMarkers)){
  ClusterMarkers[[i]]=cbind(ClusterMarkers[[i]],i)
  ClusterMarkers[[i]]=cbind(ClusterMarkers[[i]],rownames(ClusterMarkers[[i]]))
}

ClusterMarkers=do.call(rbind,ClusterMarkers)
#new colnames
colnames(ClusterMarkers)[7]<-"Genname"
colnames(ClusterMarkers)[6]<-"Cluster"

#dataset without bootstrap=GDS
GDS<-as.data.table(ClusterMarkers)

#Filting for p<0.05
GDS<-GDS[p_val<0.05]

# shows how many times a gen apppear in GDS
GDS<-GDS[,GDSGenanzahl :=.N,by=.(Genname)]
#"." equals list
#.N shows sum of counts. sort for genname and show how many times gene appears.
#"by" in 172: which group shoul be sorted by?

#keep genes that appear once (for all Clusters).
GDS<-GDS[GDSGenanzahl==1]
#Here we want only the unique markergenes.


################################
####compare with bootstrap ###
################################

#import bootstrap-data (6.1)
wd = "A:/RFILES/bootstrap/22062022"
setwd(wd)

rawfiles = list.files(wd)
rawfiles#show wchich index data has for line 178
Bootraw= read.csv(rawfiles[1], row.names = 1, sep =";")
#change colnames to same colnames as GDS 
colnames(Bootraw)[8]<-"Genname"
colnames(Bootraw)[7]<-"Cluster"
Boot<-as.data.table(Bootraw)

##################################################################################
#two Data.Tables (GDS und Boot). compare both with each other.#
##################################################################################
#show me whcih gene in GDS also appears in Boot and wchich Cluster does it have in GDS
#write Cluster-number from GDS as new col in boot

Boot<-Boot[GDS,RefCluster :=i.Cluster,on=.(Genname)]
#apply Cluster-numbers from GDS to bootstrap
#map both against genenames
#NA, if genes are not found in GDS but in Boot
#add new col RefCluster to Boot. Shows which clusternumber gene had in GDS

##########################################################################################################

#number of genes per cluster fpr every bootstrap-run
ClusterCounts<-Boot[,Clustercounts :=.N,by=.(i,Cluster,RefCluster)]
#group for Bootstrap-run, Clusternumber and number of counts
# col Clustercounts shows how many times something was found
# shows how many times a cluster with the same reflcuster was found for every Run
#its not for each individual gene. Instead its for every gene with the same Cluster and Refcluster


#remove NA
ClusterCounts<-na.omit(ClusterCounts)
#delete NAs otherwise NAs will make up the majority of results.

#show if Cluster equals RefCluster
ClusterCounts$ClusterMatch<-ClusterCounts$Cluster==ClusterCounts$RefCluster
#how many times does a gene appear in Boot
ClusterCounts<-ClusterCounts[,Genanzahl :=.N,by=.(Genname)]
#show how many genes there are for each cluster per bootstrap-run
ClusterCounts<-ClusterCounts[,CountsfürCluster :=.N,by=.(Cluster,i)]
#example Bootstraprun 1, Cluster =4, CountsfürCluster=8 means, that there are 8 Markergenes in cluster 4 of Bootstrap-Run 1
# show how many genes per cluster for all Runs
ClusterCounts<-ClusterCounts[,FeaturesproCluster :=.N,by=.(Cluster)]
#Bootanzahl: show how many genes each run has
ClusterCounts<-ClusterCounts[,BootAnzahl :=.N,by=.(i)]

#show how many Genes are in the whole bootstrap-dataset
ClusterCounts<-ClusterCounts[,Featuresanzahl :=.N]

#give col new names.
colnames(ClusterCounts)[6]<-"BootstrapRun"
colnames(Boot)[6]<-"BootstrapRun"

###################################################################################
#####################Matching######################################################
##################################################################################

#show proportion of genes with the same Cluster, Refcluster and Bootstrap-run-number compared to number of genes in the same cluster and run
ClusterCounts<-ClusterCounts[,Proportion :=(Clustercounts/CountsfürCluster)*100]

#show how many clusters and refclusters are the same
ClusterCounts<-ClusterCounts[,MatchproBootstrapCluster :=.N,by=.(BootstrapRun,ClusterMatch)]
ClusterCounts<-ClusterCounts[,ProzentanrichtigenClusternproBootstrap :=(MatchproBootstrapCluster/BootAnzahl)*100]

ClusterCounts<-ClusterCounts[,MatchproCluster :=.N,by=.(Cluster,ClusterMatch)]
ClusterCounts<-ClusterCounts[,ProzentanrichtigenClustern :=(MatchproCluster/FeaturesproCluster)*100]

# show how many genes were mapped correctly (on gene-vele)
#example: jup appears 176 times. 93 times it was assigned correclty and 83 times incorrectly
ClusterCounts<-ClusterCounts[,MatchGenEbene :=.N,by=.(Genname,ClusterMatch)]
#show percentage of times a gen was assigned to the correct refcluster.
ClusterCounts<-ClusterCounts[,ProzentRichtigeGenMatches :=(MatchGenEbene/Genanzahl)*100]

ProzentRichtigeCluster<-ClusterCounts[ClusterCounts[,.I[unique(ProzentanrichtigenClustern)],by =.(Cluster,ClusterMatch)]$V1]

##############################################################################
####################visualisation############################################
##############################################################################

########################################
######## 1  ############################
########################################

df = ProzentRichtigeCluster[,c(7,22,11)]

#add a dummy row with 0% TRUE matches for cluster 8
dummyRow = data.frame(8, 0, TRUE)
colnames(dummyRow) = colnames(df)
df = rbind(df, dummyRow)
df$ProzentanrichtigenClustern = round(df$ProzentanrichtigenClustern)

ggplot(data = df,
       aes(x = Cluster, y = ProzentanrichtigenClustern, fill = ClusterMatch))+
  geom_bar(stat = "identity",
           position = position_dodge()) +
  geom_text(aes(label = ifelse(ProzentanrichtigenClustern == 0, "",ProzentanrichtigenClustern)),
            vjust = 1.6,
            color = "white",
            position = position_dodge(0.9),
            size = 3.5) +
  labs(x = "Cluster",
       y = " %-Anteil Clusterzuweisung",
       fill = "Zuweisung") + 
  scale_x_continuous(breaks=1:8,
                     labels=as.character(1:8)) +
  scale_fill_discrete(labels = c("Falsch", "Richtig")) +
  # theme_minimal() + 
  theme(text = element_text(size = 10))


#########################################
######### 2 #############################
#########################################
#show how many genes were assigned correctl yor incorrectly. barplot
ggplot(ClusterCounts,aes(ClusterMatch,WahrFalsch,fill=ClusterMatch))+
  geom_bar(stat="identity",position=position_dodge())+geom_text(aes(label=round(WahrFalsch)), vjust=1.6, color="white",
                                                                position = position_dodge(0.9), size=3.5)+
  labs(x = "Richtige oder falsche Zuweisung über alle Bootstrap-Durchläufe",
       y = " %-Anteil richtiger Zuweisungen",
       fill = "Zuweisung")
 

#########################################
########## 3 ############################
#########################################

#show how many genes appear per bootstrap-run
T1 = table(ClusterCounts$BootstrapRun, ClusterCounts$Genname)

T2 = (T1>0)

hist(apply(T2, 2, sum), xlab="Anzahl Bootstrap-Runs", ylab="Anzahl Gene", main="", xlim=c(0, 100), cex.lab=1.5, cex.axis=1.5)

box()

#########################################
########## 4 ############################
#########################################

#Boxplot: percantage that a gene is assigned correctly or incorrectly on y-axis
ProzentRichtigeGenMatches<-ClusterCounts[ClusterCounts[,.I[which.max(ProzentRichtigeGenMatches)],by =.(ClusterMatch,Genname)]$V1]
ProzentRichtigeGenMatches<-ProzentRichtigeGenMatches[,Median :=median(ProzentRichtigeGenMatches),by=.(ClusterMatch)]
ProzentRichtigeGenMatches<-ProzentRichtigeGenMatches[,max :=max(ProzentRichtigeGenMatches),by=.(ClusterMatch)]
ProzentRichtigeGenMatches<-ProzentRichtigeGenMatches[,min :=min(ProzentRichtigeGenMatches),by=.(ClusterMatch)]

ggplot(ProzentRichtigeGenMatches,aes(y=ProzentRichtigeGenMatches,group=ClusterMatch))+geom_boxplot(fill=c("pink","darkred"),aes(x=ClusterMatch))+
  labs(x="Boxplot für Markergene mit falschen Zuweisungen und Boxplot für Markergene mit richtigen Zuweisungen",y="Reproduzierbarkeit in %",title="Prozentsatz an Markergenen im Bootstrap mit denselben Clustern wie im Gesamtdatensatz")


######################################
########## 5##########################
######################################


df<-ClusterMarkers

df<- subset(df,p_val<0.05)

C1<-df[df$Cluster == 1,7]

C2<-df[df$Cluster == 2,7]

C3<-df[df$Cluster == 3,7]

C4<-df[df$Cluster == 4,7]

C5<-df[df$Cluster == 5,7]

C6<-df[df$Cluster == 6,7]

C7<-df[df$Cluster == 7,7]

C.all = list(C1, C2, C3, C4, C5, C6, C7)


Genenames = table(df$Genname)

length(which(Genenames == 1)) # -> 214 Gene names show up

ClusterNames = paste0("Cluster ", 1:7, "\n", "(n = ", sapply(C.all, function(x) length(x)), ")")

venn(C.all, zcolor="style", snames = ClusterNames, box = FALSE, ilcs = 1, sncs = 1)
