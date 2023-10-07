###############################################################################
######################csv import###############################################
###############################################################################

filename = file.choose()
boot = read.csv(filename,sep=";",row.names=NULL)

##############Anzahl an Genen und unique Genen vergleichen#####################
nrow(boot)

length(unique(boot$rownames.iteration..i....j...))

###############################################################################

#Tabelel mit Genanzahl in allen ITerationen.Absteigende Reihenfolge
sort(table(boot$rownames.iteration..i....j...),decreasing =TRUE )

#neues Objekt "boot.p_adj": Nur Gene, die p_adj<0,05 sind.
#filtere alles raus, was P-adj größer 0,05 ist.
boot.p_adj=boot[boot$p_val_adj<=0.05,]

###############################################################################

#Tabelle auf Gene(V2)und Iterationen (V1) reduzieren
abc=as.data.frame(cbind(boot.p_adj$i,boot.p_adj$rownames.iteration..i....j...))
#angeben, wie oft Gene signifikant aufgetaucht sind
sig.genes=sort(table(unique(abc)[,2]),decreasing = TRUE)
#Darstellung verbessern
sig.genes=cbind(names(sig.genes),as.numeric(sig.genes))

###############################################################################
#############Matrix mit signifikanten Markergenen exportieren##################
###############################################################################
sig.genes=as.data.frame(sig.genes)
sig.genes$V2=as.numeric(sig.genes$V2)

write.table(sig.genes,"sigGene.csv")
