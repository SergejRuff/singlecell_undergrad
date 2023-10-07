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
library(ggplot2)#Für den Fall,dass ich Graphen mit ggplot mache.
library(patchwork)
library(venn)
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
####Gesamtdatensatz filtern##########
#####################################

ClusterMarkers=as.vector(active_clusters,mode="list")

for (i in 1:length(ClusterMarkers)){
  ClusterMarkers[[i]]=cbind(ClusterMarkers[[i]],i)
  ClusterMarkers[[i]]=cbind(ClusterMarkers[[i]],rownames(ClusterMarkers[[i]]))
}

ClusterMarkers=do.call(rbind,ClusterMarkers)
#neue Spaltennamen
colnames(ClusterMarkers)[7]<-"Genname"
colnames(ClusterMarkers)[6]<-"Cluster"

#Gesamtdatensatz=GDS
GDS<-as.data.table(ClusterMarkers)

#Filter für p<0.05
GDS<-GDS[p_val<0.05]

#wie oft taucht Gen in GDS auf?
GDS<-GDS[,GDSGenanzahl :=.N,by=.(Genname)]
#"." in 172 ist Abkürzung für Liste
#.N gibt die Anzahl (Counts) raus. Sortiere hier nach Genname und Gebe Anzahl raus.
#"by" in 172: nach welchen Gruppen soll sortiert werden?

#Behalte nur die Gene, die 1 einziges Mal auftauchen (über alle Cluster).
GDS<-GDS[GDSGenanzahl==1]
#Soll man es nur bei GDS machen oder auch alle Gene in 100 Bootstraps auch?
#Hier will man nur die unique MArkergene haben. Macht es Sinn Bootstrap zu filtern?
#Nur im GDS alle Gene rausschmeißen,die nicht unique sind oder auch Bootstrap?

################################
####Mit Boot_Strap abgleichen###
################################

#Bootstrap aufrufen
wd = "A:/RFILES/bootstrap/22062022"
setwd(wd)

rawfiles = list.files(wd)
rawfiles#überprüfe welchen Index die Datei hat für Zeile 123
Bootraw= read.csv(rawfiles[1], row.names = 1, sep =";")
#Spalte unbennen für Zeile 204, damit GDS und Boot gleichen Spaltennamen haben.
colnames(Bootraw)[8]<-"Genname"
colnames(Bootraw)[7]<-"Cluster"
Boot<-as.data.table(Bootraw)
#Boot<-Boot[p_val<0.05]

##################################################################################
#Habe jetzt zwei Data.Tables (GDS und Row). Möchte beide miteinander vergleichen.#
##################################################################################
#Will sehen, welches Gen in meinem GDS finde ich auch in Boot wieder und in welchem Cluster ist es im GDS.
#Schreibe die Clusternummer als neue Spalte in Boot.

Boot<-Boot[GDS,RefCluster :=i.Cluster,on=.(Genname)]
#Übergebe hier Bootstraps meinen Gesamtdatensatz (GDS)
#Mappe beide gegeneinander basierend auf Gennamen.
#Es sucht also Gennamen raus und mappt gegen Gennamen.
#NA, wenn Gene nicht in GDS vorzufinden sind, aber in Boot schon.
#Fügt neue Spalte RefCluster hinzu. Zeigt an, welchem Cluster es im GDS enspricht.

##########################################################################################################


#Anzahl der Gene pro Cluster für jeden Bootstraprun
ClusterCounts<-Boot[,Clustercounts :=.N,by=.(i,Cluster,RefCluster)]
#Gruppiere nach Bootstrap,Clusterzahl und nach Anzahl an gefundenen Counts.
#Clustercounts Spalte sagt, wie oft etwas vorkommt.
#Zeigt mir an, wie oft ein Cluster gefunden wurde für jeden Durchlauf verglichen mit selben RefCluster und Cluster.
#Es ist nicht pro Gen, sondern für alle Gene im selben Durchlauf, mit dem selben Cluster und dem selben RefCluster

#entferne NA
ClusterCounts<-na.omit(ClusterCounts)
#muss NAs entfernen, da ansonsten nur noch NAs übrig bleiben.

#enstpricht Cluster RefCluster?
ClusterCounts$ClusterMatch<-ClusterCounts$Cluster==ClusterCounts$RefCluster


#Wie oft kommt ein Gen vor
ClusterCounts<-ClusterCounts[,Genanzahl :=.N,by=.(Genname)]

#Wie viele sind im jeweiligen Cluster pro bootstrap-Durchlauf vorzufinden
ClusterCounts<-ClusterCounts[,CountsfürCluster :=.N,by=.(Cluster,i)]
#Bsp. Bootstraprun 1, Cluster =4, CountsfürCluster=8 bedeutet, dass es im Bootstrap-Run 1, Cluster 4 8 Markergene gibt

#Wie viele Markergene pro Cluster (über alle B-Runs)
ClusterCounts<-ClusterCounts[,FeaturesproCluster :=.N,by=.(Cluster)]

#Wie viele Markergene pro BootstrapRun
#Bootanzahl: Wie viele Rows sind pro Bootstrap enthalten
ClusterCounts<-ClusterCounts[,BootAnzahl :=.N,by=.(i)]

#Gesamte Anzahl an Features über alle BootstrapDurchläufe
ClusterCounts<-ClusterCounts[,Featuresanzahl :=.N]

#Gebe i  einen neuen Namen.
colnames(ClusterCounts)[6]<-"BootstrapRun"
colnames(Boot)[6]<-"BootstrapRun"

###################################################################################
#####################Matching######################################################
##################################################################################


#P-values aus GDS auf Bootstrap-Matrix übertragen
ClusterCounts<-ClusterCounts[GDS,RefPValue :=i.p_val,on=.(Genname)]

# Was ist die Proportion der Anzahl der Markergene mit demselben Cluster,refCluster und Bootstraprun im Verhältnis zur Anzahl an Markergenen in dem Cluster und Bootstrap 
ClusterCounts<-ClusterCounts[,Proportion :=(Clustercounts/CountsfürCluster)*100]

#Wie viele cluster und RefCluster stimmen pro Bootstrap überein?
ClusterCounts<-ClusterCounts[,MatchproBootstrapCluster :=.N,by=.(BootstrapRun,ClusterMatch)]
ClusterCounts<-ClusterCounts[,ProzentanrichtigenClusternproBootstrap :=(MatchproBootstrapCluster/BootAnzahl)*100]

ClusterCounts<-ClusterCounts[,MatchproCluster :=.N,by=.(Cluster,ClusterMatch)]
ClusterCounts<-ClusterCounts[,ProzentanrichtigenClustern :=(MatchproCluster/FeaturesproCluster)*100]

#Auf Ebene des einzelnen Gens. Wie viel vom gen wurde richtig gemappt udn wie viele falsch
#Bsp.Gen Jup kommt 176 mal vor und er wurde 93 richtig und 83 mal falsch gemappt
ClusterCounts<-ClusterCounts[,MatchGenEbene :=.N,by=.(Genname,ClusterMatch)]
#Wie oft in Prozent wurde das Gen falsch oder richtig gemappt.
ClusterCounts<-ClusterCounts[,ProzentRichtigeGenMatches :=(MatchGenEbene/Genanzahl)*100]

#Wie viele wahre und falsche Treffer gibt es
ClusterCounts<-ClusterCounts[,Wahr :=.N,by=.(ClusterMatch)]
#Wie viel % der Zuordnungen sind wahr oder falsch über die gesamte Datei/Matrix
ClusterCounts<-ClusterCounts[,WahrFalsch :=(Wahr/Featuresanzahl)*100]

#Median,Min,Max von Features pro Cluster
ClusterMedian<-ClusterCounts[ClusterCounts[,.I[unique(CountsfürCluster)],by =.(Cluster,BootstrapRun)]$V1]
ClusterMedian<-ClusterMedian[,Median :=median(CountsfürCluster),by=.(Cluster)]
ClusterMedian<-ClusterMedian[,max :=max(CountsfürCluster),by=.(Cluster)]
ClusterMedian<-ClusterMedian[,min :=min(CountsfürCluster),by=.(Cluster)]

ProzentRichtigeCluster<-ClusterCounts[ClusterCounts[,.I[unique(ProzentanrichtigenClustern)],by =.(Cluster,ClusterMatch)]$V1]
###########################################################################
######Export der Daten#####################################################
###########################################################################

#Ordner, wo gespeichert werden soll
setwd("A:/RFILES/Neuer Ansatz/Gesamtdatensatz")

#Speichere gesamtdatensatz ClusterMarkers
write.table(ClusterMarkers,"Gesamtdatensatz.csv",quote=FALSE,sep=";")

#speichere Gesamtdatensatz GDS
write.table(GDS,"GDS.csv",quote=FALSE,sep=";")

#speichere Boot
write.table(Boot,"Boot.csv",quote=FALSE,sep=";")

#speichere ClusterCounts
write.table(ClusterCounts,"ClusterCounts.csv",quote=FALSE,sep=";")

#speichere ClusterDaten
write.table(ClusterDaten,"ClusterDaten.csv",quote=FALSE,sep=";")

#speichere MaxCounts
write.table(MaxCounts,"MaxCounts.csv",quote=FALSE,sep=";")


##############################################################################
####################Visualisuerung############################################
##############################################################################

########################################
#Abbildung 1############################
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
#Abbildung 2#############################
#########################################
#Wie viele Features sind richtig oder falsch gemappt wurden über alle verfügbaren Features.
ggplot(ClusterCounts,aes(ClusterMatch,WahrFalsch,fill=ClusterMatch))+
  geom_bar(stat="identity",position=position_dodge())+geom_text(aes(label=round(WahrFalsch)), vjust=1.6, color="white",
                                                                position = position_dodge(0.9), size=3.5)+
  labs(x = "Richtige oder falsche Zuweisung über alle Bootstrap-Durchläufe",
       y = " %-Anteil richtiger Zuweisungen",
       fill = "Zuweisung")
#Problem: Warum sind die Werte so schlecht: Liegt es daran, dass ich einen Fehlergemacht habe oder liegt es an Single Cell und Bootstrapping  

#########################################
#Abbildung 3#############################
#########################################

#Wie viele Markergene tauchen pro Bootstrap-Run auf
T1 = table(ClusterCounts$BootstrapRun, ClusterCounts$Genname)

T2 = (T1>0)

hist(apply(T2, 2, sum), xlab="Anzahl Bootstrap-Runs", ylab="Anzahl Gene", main="", xlim=c(0, 100), cex.lab=1.5, cex.axis=1.5)

box()

#########################################
#Abbildung 4#############################
#########################################

#Boxplot: Wahrscheinlichkeit, dass ein Gen richtig oder falsch gemappt wird auf Y.
ProzentRichtigeGenMatches<-ClusterCounts[ClusterCounts[,.I[which.max(ProzentRichtigeGenMatches)],by =.(ClusterMatch,Genname)]$V1]
ProzentRichtigeGenMatches<-ProzentRichtigeGenMatches[,Median :=median(ProzentRichtigeGenMatches),by=.(ClusterMatch)]
ProzentRichtigeGenMatches<-ProzentRichtigeGenMatches[,max :=max(ProzentRichtigeGenMatches),by=.(ClusterMatch)]
ProzentRichtigeGenMatches<-ProzentRichtigeGenMatches[,min :=min(ProzentRichtigeGenMatches),by=.(ClusterMatch)]
#ProzentRichtigeGenMatches<-na.omit(ProzentRichtigeGenMatches)
ggplot(ProzentRichtigeGenMatches,aes(y=ProzentRichtigeGenMatches,group=ClusterMatch))+geom_boxplot(fill=c("pink","darkred"),aes(x=ClusterMatch))+
  labs(x="Boxplot für Markergene mit falschen Zuweisungen und Boxplot für Markergene mit richtigen Zuweisungen",y="Reproduzierbarkeit in %",title="Prozentsatz an Markergenen im Bootstrap mit denselben Clustern wie im Gesamtdatensatz")


######################################
#Abbildung 5##########################
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



#################################################


sigGenes<-GDS[order(p_val), head(.SD, 5), by = Cluster]

t3<-apply(T2, 2, sum)
t4<-as.data.frame(t3)
sigGene<-t4
sigGene$Genname<-rownames(sigGene)
sigGene<-as.data.table(sigGene)
sigGenes<-sigGenes[sigGene,InwievielenBootRuns :=i.t3,on=.(Genname)]

# 2 der signifikantesten Gene in GDS kommen in weniger als 20 runs vor

write.table(sigGenes,"sigGenesimGDS.csv",quote=FALSE,sep=";")
