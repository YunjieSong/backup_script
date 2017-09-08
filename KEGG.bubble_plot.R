args = commandArgs(trailingOnly=T) 

# args="dmr_KEGG.xls"

if(length(args) != 1){
  
  cat(
    "Usage: 
    Rscript bubble.r dmr_KEGG.xls
    ")
  options("show.error.messages" = F) 
  stop()
  
}

file = normalizePath(args[1])

library(ggplot2)
library(reshape2)

input = read.table(file,header=T,sep="\t")

if(length(input)>50) pathway<-input[1:50,] else  pathway<-input

a<-pathway$Size
b<-pathway$Count

richFactor=b/a

pr = ggplot(pathway,aes(richFactor,Term)) + geom_point(aes(size=Count,color=-1*log10(Pvalue))) + scale_colour_gradient(low="green",high="red") + labs(color=expression(-log[10](Pvalue)),size="Gene number",x="Rich factor",y="Pathway name",title="Statistics of Pathway Enrichment")

pr + theme_bw()

ggsave("PathWayEnrichment.pdf")   

