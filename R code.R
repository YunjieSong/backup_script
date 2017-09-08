source("http://bioconductor.org/biocLite.R")
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
biocLite("KEGG.db")

#### PCA + ggplot

pca.res=prcomp(t(maqc_dat))
pc.var=pca.res$sdev^2
pc.per=round(pc.var/sum(pc.var)*100,1)

pro1 <- paste("PC1 (",pc.per[1],"%)",sep="")
pro2 <- paste("PC2 (",pc.per[2],"%)",sep="")

pc=as.data.frame(pca.res$x)

pc$Sample=substr(rownames(pc),7,7)
pc$Site=substr(rownames(pc),5,5)
pc$Platform=substr(rownames(pc),1,3)

library(ggplot2)
PCAplot=ggplot(pc,aes(PC1,PC2))+
geom_point(aes(color=Platform,shape=Sample,size=Site))+
labs(x=pro1,y=pro2,title="PCA")

png("20160328.maqc.1og2.intensity.PCA.png",width=7*ppi,height=6*ppi,res=ppi)
PCAplot
dev.off()
####




## extract string between two words:
library(qdapRegex)
rm_between(x, "PRODUCT", "OKAY", extract=TRUE)



######  https://support.bioconductor.org/p/74303/
### If you are a user of biomaRt (a part of the Bioconductor library) change the 
### host from 'www.biomart.org' to 'www.ensembl.org'. Similar procedure can be used for all the other community servers.
library(biomaRt)
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host = "jul2015.archive.ensembl.org")

listMarts(host="www.ensembl.org")
# or
mouse <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "may2012.archive.ensembl.org")
# or
human <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "grch37.ensembl.org")
# mart=useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "grch37.ensembl.org")

##### hclust with color bar  1):
library(WGCNA)
sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE);
plotDendroAndColors(sampleTree2, traitColors,groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")

?labeledHeatmap


##### hclust with color bar  2):
library(dendextend)
plot(hclust,xlab="",sub="",hang=-1)
colored_bars(colors=dataAnnot,dend=hclust,sort_by_labels_order=T,cex.rowLabels=0.8)




###########
mycol=RColorBrewer::brewer.pal(n=7,name="Set1")


########## topGO : GO enrichmeant
library("topGO")
tgd <- new( "topGOdata", ontology="BP", allGenes = allg, nodeSize=5,annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl" )
## run tests
resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
tab<- GenTable( tgd, Fisher.elim = resultTopGO.elim,Fisher.classic = resultTopGO.classic,orderBy = "Fisher.classic" , topNodes = 200)


### gene length
library(geneLenDataBase)
genes <- c("ENSG00000124208", "ENSG00000182463", "ENSG00000124201", "ENSG00000124205", "ENSG00000124207")
getlength(genes,'hg19','ensGene')

### get GO
go=getgo(names(genes),"hg19","ensGene")

### get GOTERM

goterm=GOTERM[[go]]


### pvca
library(pvca)
mx <- all.res
colnames(mx) <- 1:ncol(all.res)
phe <- data.frame(Samples=1:ncol(all.res),tissue=factor(spe),species=factor(tis))
rownames(phe) <- 1:ncol(all.res)
library(Biobase)
input <- ExpressionSet(mx , phenoData=AnnotatedDataFrame(data=phe))
batch.factors <- c("tissue","species")
pvcaObj <- pvcaBatchAssess (input,batch.factors,threshold=0.6)


## t.test error solution
pvalue <- function(x,mm,nn) {
    obj<-try(t.test(x[1:mm],x[(mm+1):nn],var.equal = TRUE), silent=TRUE)
   if (is(obj, "try-error")) return(NA) else return(obj$p.value)
   }


## R layout
layout(matrix(c(1,2,1,2), 2, 2, byrow = TRUE), 
       widths=c(3,1))
par(mar=c(5,4,4,1))
plot(1)

par(mar=c(5,0,4,1))
plot(1,axes=F,type="n")

legend("topleft",legend="A",bty="n")



## add colors to hclust
### clustering
dist=dist(t(data))
hclust=hclust(dist)
col=colors()[1:100]

library(dendextend)
dend <- as.dendrogram(hclust)
labels_colors(dend) <- col[order.dendrogram(dend)]
plot(dend)


##使用BiocParallel来进行并行运算
library("BiocParallel")
library(BatchJobs)
cluster.functions <- makeClusterFunctionsMulticore()
# Setting for worker localhost: ncpus=8
bpParams <- BatchJobsParam(cluster.functions=cluster.functions)
system.time(r2 <- bplapply(buf[1:10000], parser, BPPARAM=bpParams))



## solve x according to Y
m=lm(x~y)
m$coefficients%*%matrix(c(1,50),ncol=1)

###    Principal Variance Component Analysis 
library(pvca)
?pvcaBatchAssess

library(VennDiagram)
venn.diagram(list(GSE12093=re.12093,GSE2990=re.2990,GSE6532=re.6532),"tamoxi_genes.png",fill=c("red","blue","yellow"))


##  compute permutation adjusted p-values form multiple testing
## library(multtest)
mt.maxT()

# FTP 数据下载
lftp -c "pget -n 10 ftp://ftp-private.ncbi.nih.gov/SEQC/tmp/Leming/D0C8LACXX_748_8P2_TGACCA_L005_R2_001.fastq.gz"
lftp -c "pget -n 2 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn/GSE25065/matrix/GSE25065_series_matrix.txt.gz"
# from GEO
lftp -c "pget -n 10 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19615/suppl/GSE19615_RAW.tar"
# or
library(GEOquery)
tmp=getGEOSuppFiles("GSE12093")


# find the position of the max value of a matrix
which(data==max(data),arr.ind=T)  # "data" is a numeric matrix 

# return the position of the upptriange position of a matrix
res=c()
for(i in 2:nrow(data))
{
for(j in 1:(i-1))
{
tmp=paste(j,i,sep="_")
res=c(res,tmp)
}
}
# data is matrix

# reduced HCA label size
par(cex=0.3, mar=c(5, 8, 4, 1))
plot(hc, xlab="", ylab="", main="", sub="", axes=FALSE)
par(cex=1)
title(xlab="xlab", ylab="ylab", main="main")
axis(2)


#正态性检验
NormTest<-function(data){  
library(fBasics)  
library(nortest)  
udata<-unique(data)  
result<-list()  
result$D<-dagoTest(data)  
result$jB<-jarqueberaTest(data)  
result$SW<-shapiroTest(data)  
result$lillie<-lillie.test(data)  
result$ad<-ad.test(data)  
result$cvm<-cvm.test(data)  
result$sf<-sf.test(data)  
return(result)  
}

###
model=lm(weight ~ group)
library(car) 
outlierTest(model)
influence.measures(model)

## upper.tri()

##read GEO matrix
exp=read.table("GSE48391-GPL570_series_matrix.txt",skip=67,sep="\t",fill=NA,header=T,row.names=1)


##microarray array probe annotation:
biocLite("hgu133plus2.db")
library(hgu133plus2.db)
library(annotate)
gene.symbols <- getSYMBOL(probeset.list$ID, "hgu133plus2")


## snp.plotter
library("snp.plotter")


# remove missing values 
# data = uis_samll
tiny_uis<-uis_small[apply(uis_small,1,function(x)!any(is.na(x))),]

====================================
CV:
var.coeff
function (x, square = FALSE) 
{
n <- length(x)
V <- sqrt((n - 1) * var(x)/n)/mean(x)
if (square) 
V <- V^2
V
}

coev(变异系数)

coe.v=function(x){sd(x)/mean(x)}
=====================================
Secondary Y axis title(第二坐标)

par(new = TRUE)

plot(z[,2], type = "l", ann = FALSE, yaxt = "n", col = "blue", lty = 2, lwd 
= 2)

axis(4)

legend(x = "topleft",  bty = "n",  lty = c(1,2),  lwd = c(2,2),
        col = c("black", "blue"),  legend = paste(colnames(z), c("(left 
y-axis)", "(right y-axis)")))

mtext(colnames(z)[2], 4, line = 3, col = "blue")

-----------------------------------------------------------------
time=c(2,4,10)
y1=c(9.2,7.667,4.667)
y2=c(4.396,19.082,73.984)
plot(time,y1,type='l',col='red',ylab="Relative levels", xlab="Age of Mice",ylim=c(0,20)) 
points(time,y1,col="red", pch=19)

op <-par(new=T) 
plot(time,y2,type='l', col='green',axes=F,xlab="",ylab="", ylim=c(0,80))
points(time,y2,col='green',pch=24) 
axis(4)
----------------------------------------------------------
#################(secondary y axis)
x <- 1:5
y1 <- rnorm(5)
y2 <- rnorm(5,20)
par(mar=c(5,4,4,5)+.1)
plot(x,y1,type="l",col="red")
par(new=TRUE)
plot(x, y2,,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="")
axis(4)
mtext("y2",side=4,line=3)
legend("topleft",col=c("red","blue"),lty=1,legend=c("y1","y2"))

##############################################read table from website
library(XML)
tables=readHTMLTable("http://www.fda.gov/drugs/scienceresearch/researchareas/pharmacogenetics/ucm083378.htm")
tables
library(xlsx)
write.xlsx(tables,"tables.xlsx")

#########
划分窗口及设置页边距
windows(width=20,height=10)
par(1:2,mfrow=c(1,2)，mar=c(6,5,5,4))
######################################################
设置坐标轴
plot(1:4, rnorm(4), axes = FALSE)
axis(1, 1:4, LETTERS[1:4])
axis(2)
====================================
改变坐x标轴的标记
plot(x,y,xaxt="n")
axis(side=1,at=c(....),labels=c(...))
或
x=c(1:4)
y=c(5:8)
boxplot(x,y,xaxt="n")
axis(side=1,at=c(1,2),labels=c("controlA","controlB"))
====================================
vennDiagram(维恩图)
draw.pairwise.venn (4794,2007,1704,fill=c("red","blue"),category = c("Anova_FC", "T_FC"),cex = 2)


1 library(plotrix)  

2 plot(0:10,seq(0,10,length=11),type='n',axes=F,xlab='',ylab='') draw.circle(2,5,2,col=rgb(154/255,0/255,205/255,0.6))  

3 draw.circle(4,5,2,col=rgb(21/255,3/255,252/255,0.6)) text(1,5,labels='10.12%',col='white',font=2) text(5,5,labels='40.38%',col='white',font=2)  

4 text(3,5,labels='49.5%',col='white',font=2)  

5 legend(6.2,5,pch=15, xjust=0,yjust=0.5,bty='n',cex=1.3,col=c(rgb(154/255,0/255,205/255),rgb(74/255,2/255,233/255),rgb(21/255,3/255,252/255) ), legend=c('sample 1 uniq','sample& sample 2','sample 2 uniq'))  

6 text(3.5,7.5,labels='Venn chart for uniq_sRNAs',font=2,cex=1.5) 






text(locator(1),"A",srt(文字方向)，xpd=T (#文字区域))


=========================================两两比对pvalue


pvalue=function(x,mm,nn) {
    obj<-try(t.test(x[1:mm],x[(mm+1):nn],var.equal = TRUE), silent=TRUE)
   if (is(obj, "try-error")) return(c(1)) else return(obj$p.value)
   }


diff.exp=matrix(,ncol=?,nrow=?)  
k=0
AB=c()
for (i in 1:4){
for (j in (i+1):5) {
exp1=EXP[,(3*j-2):(3*j)]
exp2=EXP[,(3*i-2):(3*i)]
exp=cbind(exp1,exp2)
mm=ncol(exp1)
nn=ncol(exp)
k=k+1
diff.exp[,k]=apply(exp,1,function(x)pvalue(x=x,mm=mm,nn=nn))
}
}




#####################原代码：
ttest.mvsn<-matrix(0,ncol=153,nrow=dim(raw)[1])
nam<-c()
n=0
for (i in 1:17){
	for (j in (i+1):18){
o<-paste(utis[i],utis[j],sep="_")
exp1<-raw[,grep(paste(utis[i],"_",sep=""),colnames(raw))]
exp2<-raw[,grep(paste(utis[j],"_",sep=""),colnames(raw))]
mm<-ncol(exp1)
exp<-cbind(exp1,exp2)
nn=ncol(exp)
n<-n+1
ttest.mvsn[,n]<-apply(exp,1,function(x)pvalue(x=x,mm=mm,nn=nn))
nam<-c(nam,o)
print (o)
}
}



或者 （张立）

AB=matrix(,ncol=10,nrow=dim(EXP)[1])
my.ttest=function(data)
{
sample1=data[1:3]
sample2=data[4:6]
p.value=t.test(sample1,sample2)$p.value
return(p.value)
}
k=0
for (i in 1:4){
for (j in (i+1):5) {
k=k+1
data1=EXP[,(i*3-2):(3*i)]
data2=EXP[,(j*3-2):(3*j)]
data=cbind(data1,data2)
print(head(data))
AB[,k]=matrix(apply(data,1,my.ttest),ncol=1)
}
}

==================================

k=0
AB=c()
for (i in 1:4){
for (j in (i+1):5) {
k=k+1
AB[k]=t.test(EXP[1,(3*j-2):(3*j)],EXP[1,(3*i-2):(3*i)])$p.value
}
}

t.test(EXP[1,4:6],EXP[1,1:3])$p.value
t.test(EXP[1,7:9],EXP[1,1:3])$p.value
t.test(EXP[1,10:12],EXP[1,1:3])$p.value
t.test(EXP[1,13:15],EXP[1,1:3])$p.value


####++++++++++++++++++多探针对应一个基因，返回最大值。。。（张立）


data=matrix(abs(rnorm(60,0,4)),ncol=10)
rownames(data)=c("a","a","b","c","d","d")
##data is 6*10 matrix

genename.set=unique(rownames(data))
probeset=paste("probe",c(1:dim(data)[1]),sep="")
max.data=apply(data,1,max)
unique.probes=c()
for(i in 1:length(genename.set))
{
 fold.change=max.data[rownames(data)==genename.set[i]]
 probes=probeset[rownames(data)==genename.set[i]]
 unique.probes[i]=probes[abs(fold.change)==max(abs(fold.change))]
}


position=match(unique.probes,probeset)


################(另一代码)重复基因求均值  (师姐)
data<-data.frame(tapply(exp[,1],factor.name,mean))
colnames(data)<-colnames(exp)[1]
for (i in 2:4) {
b<-data.frame(tapply(exp[,i],factor.name,mean))
colnames(b)<-colnames(exp)[i]
data<-merge(data,b,by=0)
rownames(data)<-data[[1]]
data[[1]]<-NULL}
####################################





### merge的用法



##多重条形图
barplot(x,beside=T，names.arg=1:5) 

#x为矩阵    
##segments()可用于绘制误差线。
#names.arg=   可用来设定横坐标 
 
#     或者    b=barplot(x,beside=T，names.arg=1:5) 
#             axis(1,at=b,labels="")



#行列数目不一致的矩阵合并：(张立)
   combine=function(big,small)
{
  rownames.big=rownames(big)
  rownames.small=rownames(small)
  union.row=union(rownames.big,rownames.small)
  combine=matrix(,nrow=length(union.row),ncol=ncol(big)+ncol(small))
  rownames(combine)=union.row
  combine[rownames.big,1:ncol(big)]=big
  combine[rownames.small,(ncol(big)+1):(ncol(big)+ncol(small))]=small
  
  return(combine)
}


##不同平台基因芯片比对（mouse to human,张立）

data=read.table("homologene.data",sep="\t")
data=as.matrix(data)
col1=as.numeric(data[,2])
human.rows=col1==9606
human.data=data[human.rows,c(1,4)]
mouse.rows=col1==10090
mouse.data=data[mouse.rows,c(1,4)]
human.upper=toupper(paste(human.data[,1],human.data[,2],sep="_"))
mouse.upper=toupper(paste(mouse.data[,1],mouse.data[,2],sep="_"))
overlap=intersect(human.upper,mouse.upper)
n=length(human.upper)+length(mouse.upper)-length(overlap)
unique.human=setdiff(human.upper,overlap)
unique.mouse=setdiff(mouse.upper,overlap)
result=matrix(,ncol=3,nrow=length(union(human.upper,mouse.upper)))
trans=function(s)
{
  n=length(s[[1]])
  m=length(s)
  res=matrix(,ncol=n,nrow=m)
  for(i in 1:m)
  {
    res[i,]=s[[i]]
  }
  return(res)
}
result[1:length(unique.human),1:2]=trans(strsplit(unique.human,"_"))
result[(length(unique.human)+1):(length(unique.human)+length(unique.mouse)),c(1,3)]=trans(strsplit(unique.mouse,"_"))
overlap1=trans(strsplit(overlap,"_"))
col2=overlap1[,2]
col.first=substr(col2,1,1)
col.res=tolower(substr(col2,2,100))
overlap2=paste(col.first,col.res,sep="")
overlap3=cbind(overlap1,overlap2)
result[(length(unique.human)+length(unique.mouse)+1):nrow(result),]=overlap3


merge(a,b,by.x="a",by.y="a",all=T)



## 定义函数，用均值代替缺失值
replace.missing=function(x)
{
x[is.na(x)]=mean(x[!is.na(x)])
return (x)
}

#定义函数Concordance Correlation Coefficient（CCC）
CCC=function(aa,bb)
{
meana=mean(aa)
meanb=mean(bb)
cor=cor(aa,bb)
sda=sd(aa)
sdb=sd(bb)
return(2*cor*sda*sdb/sum(sda^2,sdb^2,(meana-meanb)^2))
}


library(plotrix)
slices <- c(10, 12, 4, 16, 8) 
lbls <- c("US", "UK", "Australia", "Germany", "France")
pie3D(slices,labels=lbls,explode=0.2,theta=0.7,
   main="Pie Chart of Countries ")

#对奇数位置求和
a=function(x)
{
n=length(x)
b=c()
for(i in 1:n){
if(i%%2==1){
b=c(b,x[i])
}
}
return(sum(b))
}

##猴子当大王
houzi<-function(x){
x1<-data.frame(c(1:x))
num<-rep(c(1:3),x)
x1[,2]<-num[1:x]
f<-x1[which(!x1[,2]==c(3)),]
while(dim(f)[1]>=3){
f<-x1[which(!x1[,dim(x1)[2]]==c(3)),]
x1<-x1[nrow(x1):1, ]
x1[which(!x1[,dim(x1)[2]]==c(3)),dim(x1)[2]+1]<-num[1:(dim(f)[1])]
}
return(rownames(x1)[which(x1[,dim(x1)[2]]==1)])
}

一个例子

本例来自于"introduction to Scientific Programming and Simulatoin Using R"一书的习题。有这样一种赌博游戏，
赌客首先将两个骰子随机抛掷第一次，如果点数和出现7或11，则赢得游戏，游戏结束。如果没有出现7或11，赌客继续抛掷，
如果点数与第一次扔的点数一样，则赢得游戏，游戏结束，如果点数为7或11则输掉游戏，游戏结束。如果出现其它情况，则继续抛掷，
直到赢或者输。用R编程来计算赌客赢的概率，以决定是否应该参加这个游戏。


craps <- function() {
#returns TRUE if you win, FALSE otherwise
initial.roll <- sum(sample(1:6,2,replace=T))
if (initial.roll == 7 || initial.roll == 11) return(TRUE)
while (TRUE) {
current.roll <- sum(sample(1:6,2,replace=T))
if (current.roll == 7 || current.roll == 11) {
return(FALSE)
} else if (current.roll == initial.roll) {
return(TRUE)
}
}
}
mean(replicate(10000, craps()))


##linux 下从本地安装R包
setwd("./")     #更改工作目录至linux安装包所在的文件夹位置
install.packages("jetset_1.4.0.tar.gz", repos=NULL, type="source")    #“jetset_1.4.0.tar.gz”为安装包的名字

##换行
plot(rnorm(100), main="this is my title \non two lines")

a=cat("Identical cell lines\nCGP vs CCLE\n")

expression(atop("Histogram of "*hat(mu), Bootstrap~samples*','~Allianz))


plot(1,2,xlab=expression(Example:~~H[2]~ABS[3]~m^2~s^3),ylab='y')

######手动添加boxplot的x轴标签（有角度x标签）

par(mar=c(8,6,3,3))
at=boxplot(common,outline=F,col=c(rep(c("lightblue","pink"),each=15)),las=2,ylab=expression(-log10(IC[50])),xaxt="n",cex.lab=1.5)

labs=at$names       ##原始矩阵纵坐标标签
axis(1, labels = FALSE)

# Plot x labs at default x position
text(x =  seq_along(labs), y = par("usr")[3]-0.2, srt = 45, adj = 1,labels = labs, xpd = TRUE)

#################correlation matrix
library(lattice)
par(mar=c(6,5,5,4))
levelplot(AA,scales=list(x=list(rot=90)),xlab="",ylab="",main="Correlation Matrix")



###加快R的运行速度
doit=function(x)(x)^2+2*x
  
# 常规方法 
system.time(res<-lapply(1:500000,doit))

#并行方法 1 
library(parallel)
detectCores(logical = F)

cl=makeCluster(getOption("cl.cores",5)) #uss 5 cores
system.time(res<-parLapply(cl,1:500000,doit))
stopCluster(cl)
# cl <- makeCluster(4)  # 初始化四核心集群
# results <- parLapply(cl,1:x,myfun) # lapply的并行版本


#######并行计算 2 
mc <- getOption("mc.cores", 8)
results <- mclapply(1:x,myfun,mc.cores=mc)

stopCluster(cl) # 关闭集群





panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
 
# Plot #2: same as above, but add loess smoother in lower and correlation in upper
pairs(~Sepal.Length+Sepal.Width+Petal.Length+Petal.Width, data=iris,
      lower.panel=panel.smooth, upper.panel=panel.cor, 
      pch=20, main="Iris Scatterplot Matrix")











