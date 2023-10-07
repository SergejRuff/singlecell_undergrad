#Plot  Counts gegen Bootstrap und vergleich welche Cluster wie viele Counts haben und welche dieser Cluster
#nicht mit der Referenz übereinstimmen.
#plot1<-ggplot(ClusterDaten,aes(x=BootstrapRun,y=Clustercounts,color=factor(Cluster)))+geom_point()+theme_bw()
#plot2<-ggplot(ClusterDaten,aes(x=BootstrapRun,y=Clustercounts,color=ClusterMatch))+geom_point()+theme_bw()
#plot1+plot2+plot_layout(nrow=2)
#Problem: Zu viele Punkte. Unübersichtlich

#Bootstrap gegen p_value
#ggplot(ClusterCounts,aes(x=BootstrapRun,y=p_val,color=factor(Cluster)))+geom_point()+theme_bw()
#ggplot(ClusterCounts,aes(x=BootstrapRun,y=RefPValue,color=factor(Cluster)))+geom_point()+theme_bw()

#Ursprünglicher Plan
#plot3<-ggplot(ClusterCounts,aes(x=Genname,y=Genanzahl,color=factor(RefCluster)))+geom_point()+theme_bw()+
#  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))
#plot4<-ggplot(ClusterCounts,aes(x=Genname,y=Genanzahl,color=ClusterMatch))+geom_point()+theme_bw()
#plot3+plot4+plot_layout(nrow=2)

#ggplot(ClusterCounts,aes(Cluster,RefCluster))+
#  geom_tile(fill=ClusterMatch)

#Wie viele der Cluster im Bootstrap stimmen auch mit dem RefCluster überein
#ggplot(ClusterCounts,aes(BootstrapRun,fill=ClustervsNewCluster))+
#  geom_bar()
#Problem: Wie kann ich die Prozentzahl visualisieren und die Bootstraps klarer abgrenzen.
#ggplot(ClusterCounts,aes(newCluster,fill=ClustervsNewCluster))+
#  geom_bar()
#Template für Boxplot
#ggplot(ClusterCounts,aes(factor(Cluster),Proportion,fill=Proportion))+
#  geom_boxplot()

ggplot(ClusterCounts,aes(x=BootstrapRun,y=Wahrheitsgehalt,color=ClusterMatch))+geom_point()+theme_bw()  


ggplot(ClusterCounts,aes(x=BootstrapRun))+geom_histogram(binwidth=1,color="white",fill="steelblue")+theme_bw() 

###########################################################
ggplot(ClusterCounts,aes(BootstrapRun,fill=ClusterMatch))+
  geom_bar()

#Anzahl an richtigen und falschen Werten
ClusterCounts$Wahrheitsgehalt<-as.numeric(as.vector(df2$ClusterCounts$Wahrheitsgehalt))
ggplot(ClusterCounts,aes(BootstrapRun,Wahrheitsgehalt,fill=ClusterMatch))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()+geom_text(aes(label=round(Wahrheitsgehalt)), vjust=1.6, color="white",
                            position = position_dodge(0.9), size=3.5)

#spearman korrelationstest
ggplot(ClusterCounts,aes(x=RefPValue,y=p_val))+geom_point()+theme_bw()+coord_trans(y = "log10",x = "log10")+geom_smooth(method=lm,formula = y ~ x, se = FALSE)  

cor(ClusterCounts$p_val,ClusterCounts$RefPValue,method="spearman")       
cor.test(ClusterCounts$p_val,ClusterCounts$RefPValue,method="spearman",exact=FALSE)

#Mittelwert pro BootstrapRun
ClusterCounts<-ClusterCounts[,Mittelwert :=mean(nrow),by=.(Cluster)]
#Abweichung pro BootstrapRun
#ClusterCounts<-ClusterCounts[,SD :=sd(Mittelwert),by=.(BootstrapRun)]


#Anzeigen, wie viele Cluster gefunden werden pro Bootstrap udn wie oft ein Cluster genutzt wird.
B<-ClusterCounts[,sum:=.N,by=.(BootstrapRun)]
B<-B[,.(X2= max(Cluster)), by = BootstrapRun]
B<-B[,sum:=.N,by=.(X2)]



#Wie viele cluster und RefCluster stimmen pro Bootstrap überein?
ClusterCounts<-ClusterCounts[,MatchCluster :=.N,by=.(Cluster,ClustervsNewCluster)]
#Bootanzahl: Wie viele Rows sind pro Bootstrap enthalten
ClusterCounts<-ClusterCounts[,ClusterAnzahl :=.N,by=.(Cluster)]
ClusterCounts<-ClusterCounts[,ProzentanrichtigenClustern :=(MatchCluster/ClusterAnzahl)*100]

ggplot(ClusterCounts,aes(Cluster))+
  geom_bar()


ggplot(ClusterCounts,aes(Cluster,ProzentanrichtigenClustern))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()+geom_text(aes(label=round(ProzentanrichtigenClustern)), vjust=1.6, color="white",
                            position = position_dodge(0.9), size=3.5)

######################################################################################

#Werden die Daten in ClusterDaten richtig zugeordnet?
ggplot(ClusterDaten,aes(Cluster,RichtigeZuweisunginProzent,fill=ClusterMatch))+
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+geom_text(aes(label=round(RichtigeZuweisunginProzent)), vjust=1.6, color="white",
                            position = position_dodge(0.9), size=3.5)+
  labs(x="Cluster",y="Richtige Zuordnung in %",title="Werden die meisten Markergene dem richtigen Referenz-Cluster zugeordnet?")
#Nur in 24 % aller BootstrapDurchläufe erhalten die meisten Markergene auch Cluster 2 als Referenz. 76 % aller Bootstraps enthalten 
#eine Mehrheit an Markergenen, die falsch zugeoodnet wurden in dem Cluster

ggplot(ClusterDaten,aes(BootstrapRun,RefCluster,color=factor(Cluster)))+
  geom_line()


ClusterDaten<-ClusterDaten[,p :=.N,by=.(Cluster,RefCluster)]
HäufigsterClusterinDaten<-ClusterDaten[ClusterDaten[,.I[which.max(p)],by =.(Cluster)]$V1]
ClusterDaten<-ClusterDaten[HäufigsterClusterinDaten,HäufigstesCluster :=i.RefCluster,on=.(Cluster)]
ClusterDaten$Match<-ClusterDaten$Cluster==ClusterDaten$HäufigstesCluster
ggplot(ClusterDaten,aes(Cluster,HäufigstesCluster,fill=Match))+
  geom_bar(stat="identity", position=position_dodge())

ggplot(MittelwertanCountsproCluster,aes(Cluster,Mittelwert))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()+geom_text(aes(label=round(Mittelwert)), vjust=1.6, color="white",
                            position = position_dodge(0.9), size=3.5)


#ClusterCounts<-ClusterCounts[,median :=median(ProzentRichtigeGenMatches),by=.(ClusterMatch)]
#################################################################################################################

##Mittelwert und SD von Counts pro Cluster für Tabelle in Bachelorarbeit

CountproCluster<-ClusterCounts[,.(CountproCluster= max(CountsfürCluster)), by = .(BootstrapRun,Cluster)]
#Mittelwert pro BootstrapRun
MittelwertanCountsproCluster<-CountproCluster[,Mittelwert :=mean(CountproCluster),by=.(Cluster)]
#Abweichung pro BootstrapRun
MittelwertanCountsproCluster<-MittelwertanCountsproCluster[,SD :=sd(CountproCluster),by=.(Cluster)]

ggplot(ClusterCounts,aes(BootstrapRun,ProzentanrichtigenClusternproBootstrap,color=ClusterMatch))+
  geom_line()


Cluster<-ClusterCounts
Cluster<-ClusterCounts[ClusterCounts[,.I[unique(Clustercounts)],by =.(BootstrapRun,Cluster,RefCluster)]$V1]
Cluster<-Cluster[,Mittelwert :=mean(Clustercounts),by=.(Cluster,Boostraprun)]
Cluster<-Cluster[,Proportion :=(Mittelwert/23760)*100]
ggplot(Cluster,aes(Cluster,Proportion,fill=factor(RefCluster)))+
  geom_bar(stat="identity", position=position_dodge())

#####################################

two<-ClusterCounts[ClusterCounts$RefCluster=="2",]
twotable<-table(two$Genname,two$Cluster)
write.table(twotable,"twotable.csv",quote=FALSE,sep=";")

#########################################################################

#enstpricht RefCluster NewCluster?
#ClusterCounts$RefClustervsNewCluster<-ClusterCounts$RefCluster==ClusterCounts$newCluster
#TRUE: RefCluster und newcluster sind identisch
#False: RefCluster und newcluster sind verschieden.

#enstpricht Cluster RefCluster?
#ClusterCounts$ClusterMatch<-ClusterCounts$Cluster==ClusterCounts$RefCluster
#TRUE: Cluster und Refcluster sind identisch
#False: Cluster und Refcluster sind verschieden.

#Wie viele cluster und RefCluster stimmen pro Bootstrap überein?
#ClusterCounts<-ClusterCounts[,MatchproBootstrapCluster :=.N,by=.(BootstrapRun,ClusterMatch)]
#Bootanzahl: Wie viele Rows sind pro Bootstrap enthalten
#<-ClusterCounts[,BootAnzahl :=.N,by=.(BootstrapRun)]
#ClusterCounts<-ClusterCounts[,ProzentanrichtigenClusternproBootstrap :=(MatchproBootstrapCluster/BootAnzahl)*100]

#enstpricht Cluster NewCluster?
#ClusterCounts$ClustervsNewCluster<-ClusterCounts$Cluster==ClusterCounts$newCluster
#TRUE: Cluster und newcluster sind identisch
#False: Cluster und newcluster sind verschieden.

#Wie viele Gene (Counts) gibt es pro Bootstrap
#ClusterCounts<-ClusterCounts[,CountsproBootsstrap :=sum(Clustercounts),by=.(BootstrapRun)]

#Proportion der Counts im Vergleich zu Counts pro Bootstrap
#ClusterCounts<-ClusterCounts[,Proportion :=(Clustercounts/CountsproBootsstrap)*100]
#Proportion von Genen in 1 Cluster proportional zur Anzahl aller Gene, 
#die in 1. Bootstraprun gefunden wurden. In %
#0.25% für Bootstraprun1, Cluster 1 und RefCluster 1 bedeutet, dass 0.25 % aller Gene,
#im Bootstrap 1 sowohl Cluster 1 als auch Refcluster 1 haben.

#Gib maximalen ClusterCount für jedes Cluster und jede Iteration raus.
MaxCounts<-ClusterCounts[ClusterCounts[,.I[which.max(Clustercounts)],by =.(i,Cluster)]$V1]
colnames(MaxCounts)[6]<-"BootstrapRun"
#behalte nur noch bootstrraprun, Clusternummer,Referenzcluster,Counts und Genname
ClusterDaten<-MaxCounts[,c(6,7,9,10,13)]
#Habe jetzt nur den Gennamen.Die Iteration.Die Clusternummer des Bootstraps. Welchem Cluster es im Referenz (GDS) entspricht und wie oft es gezählt wurde.

#Proportion der Counts im Vergleich zu Gesamtcounts pro Cluster
ClusterDaten<-ClusterDaten[,Proportion :=(Clustercounts/CountsfürCluster)*100]
#Proportion von ClusterCounts in 1 Cluster proportional zur Anzahl aller Gene, 
#die in 1. Bootstraprun gefunden wurden. In %
#Bsp. Bootstraprun=1,Cluster=3,RefCluster=3,Proportion=81% bedeutet, dass 81% aller Gene im 1. Bootstrap_Durchlauf im 3. cluster den RefCluster 3 zugeordnet bekommen haben.



#enstpricht Cluster RefCluster?
ClusterDaten$ClusterMatch<-ClusterDaten$Cluster==ClusterDaten$RefCluster
#TRUE: Cluster und Refcluster sidn identisch
#False: Cluster und Refcluster sind verschieden.

#Wie viele wahre und falsche Treffer gibt es
ClusterDaten<-ClusterDaten[,Wahr :=.N,by=.(ClusterMatch)]
#Wie viel % der Zuordnungen sind wahr oder falsch über die gesamte Datei/Matrix
ClusterDaten<-ClusterDaten[,Länge :=.N]
ClusterDaten<-ClusterDaten[,WahrFalsch :=(Wahr/Länge)*100]

#Wie viele Wahr/Falsch_Aussagen pro Cluster stimmen.
ClusterDaten<-ClusterDaten[,RichtigeZuweisungproCluster :=.N,by=.(ClusterMatch,Cluster)]
ClusterDaten<-ClusterDaten[,Clusteranzahl :=.N,by=Cluster]
ClusterDaten<-ClusterDaten[,RichtigeZuweisunginProzent :=(RichtigeZuweisungproCluster/Clusteranzahl)*100]

############################################
#Verteilung der features in RefCluster pro Cluster#
#Prozentangabe pro Spalte#
#BSP Man sieht, dass 43 % der Gene, die RefCluster 7 haben sollten, Cluster 6 zugeordnet werden.
table<-table(ClusterCounts$Cluster,ClusterCounts$RefCluster)
table2<-round(100*prop.table(table, margin = 2),digits=2)
write.table(table2,"table2.csv",quote=FALSE,sep=";")


#Übertrage die häufigsten RefCluster pro Bootstrap udn Cluster in clusterCounts zur Übersicht
ClusterCounts<-ClusterCounts[ClusterDaten,c('newCluster','Counts'):=.(i.RefCluster,i.Clustercounts),on=.(BootstrapRun,Cluster)]