args = commandArgs(trailingOnly=T) 

if(length(args) != 3){
  
  cat(
    "Usage: 
         Rscript  GO_bar_plot.R dmr_BP.xls dmr_MF.xls dmr_CC.xls
")
  options("show.error.messages" = F) 
  stop()
  
}

BP = normalizePath(args[1])
MF   = normalizePath(args[2])
CC  = normalizePath(args[3])

library(ggplot2)
library(reshape2)

read_go<-function(file,cat){
inputdata = read.table(file, header=TRUE, row.names=1, check.name=F, comment.char="", quote="", sep="\t", fill=T) 
 df = inputdata[,c(6,7)]
df$FDR = -log10(df$FDR)
df = df[1:10,]
df$category=cat
return(df)
}

df = rbind(read_go(BP, "Bological_Process"),
           read_go(MF, "Molecular_Function"),
           read_go(CC, "Cellular_Component"))


df.long=melt(df)

tiff(filename="categoryPlot.tiff",width=28,height=18,units="cm",compression="lzw",bg="white",res=600);
df.long$Term=factor(df.long$Term,levels=df.long$Term[30:1])
ggplot(df.long,aes(x=Term,y=value,fill=category,color))+geom_bar(stat="identity")+coord_flip()+ylab("-Log10(FDR)")
dev.off()


