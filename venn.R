install.packages("venn")
library("venn")
AllGenes = 1:2000
C1 = sample(AllGenes, 322, replace=FALSE)
C2 = sample(AllGenes, 324, replace=FALSE)
C3 = sample(AllGenes, 381, replace=FALSE)
C4 = sample(AllGenes, 102, replace=FALSE)
C5 = sample(AllGenes, 404, replace=FALSE)
C6 = sample(AllGenes, 410, replace=FALSE)
C7 = sample(AllGenes, 371, replace=FALSE)
C.all = list(C1, C2, C3, C4, C5, C6, C7)
venn(C.all, zcolor="style")

C1 = LETTERS[1:20]

C2 = LETTERS[15:24]

C.all = list(C1,C2)

venn(C.all, zcolor = "style")


