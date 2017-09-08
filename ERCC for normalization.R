library(DESeq)

cds = newCountDataSet( Just_ERCC_Bclass, group )
cds = estimateSizeFactors( cds )

library(edgeR)
my <- DGEList(counts=Not_ERCC, group=group)
my$samples$lib.size<-my$samples$lib.size/sizeFactors( cds )
my <- calcNormFactors(my)

## .... and so on as in the manual. 