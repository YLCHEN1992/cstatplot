cann=function(x){
up="Map R translation and annotation\n"
down="\nOver 2021.6.17\nneed packages: ggplot2; pheatmap; VennDiagram; scatterplot3d; WGCNA\n"
ncvlo="\ncvlo: need 3 cols \n    LOG->X values (log2);\n    TT->Y values (P-values);\n    TSS->point colors\n"
ncprp="\ncprp: need 3 cols \n    Bs->germ types;\n    Yv->Y values;\n    Group->Treatment type\n"
ncbarp="\ncbarp:  need 2 cols \n    Group->Treatment type; \n    Yv->Y values\n"
nchot="\nchot:  need n cols \n    col1:2->GenenamesID and <Class>; \n    cols...->samplenames\nnote: need rules\n"
ncven="\ncven: you just to fill by cols for items, and the max cols is five. So, It's freedom\n"
ncpca2d="\ncpca2d: need data.frame \n    free.sampleID->sample id need regular;\n    indexn->this is free indexs need dim>2\n"
ncpca3d="\ncpca3d: need data.frame \n    free.sampleID->sample id need regular;\n    indexn->this is free indexs need dim>3\n"
ncwgcna="\ncwgcna: need 1 or 2 files \n    1$geneID->freenames; 1$....->freesamplenames;\n    2[,1]->samplesnames==1$...;    2[,...]=traits\n"
ncgesa="\ncgesa: need 2col  \n    ENTREZID->ENTREZID;\n    LOG2->LOG2;\n    Note: need Values must sort from large to small!!!\n"
npacks="\nggplot2; pheatmap; VennDiagram; scatterplot3d; WGCNA;\n"
{cat(up);cat(get(deparse(substitute(x))));cat(down)}}

cbis=function(){install.packages("BiocManager");
BiocManager::install(c("AnnotationDbi", "impute","GO.db", "preprocessCore","org.Hs.eg.db","clusterProfiler","enrichplot"))
install.packages(c("matrixStats", "Hmisc","foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
install.packages(c("WGCNA", "stringr", "reshape2"))
install.packages(c("PerformanceAnalytics","ggplot2","pheatmap","VennDiagram","scatterplot3d","WGCNA"))}

cbi=function(){install.packages(c("BiocManager","PerformanceAnalytics","ggplot2","pheatmap",
"VennDiagram","scatterplot3d","WGCNA"))}

cbl=function(){library(ggplot2);library(pheatmap);library(VennDiagram);
library(scatterplot3d);library(WGCNA);library(PerformanceAnalytics)}


cpacks=c("BiocManager","stringr","reshape2","WGCNA","matrixStats", "Hmisc","foreach", 
"doParallel", "fastcluster", "dynamicTreeCut", "survival","AnnotationDbi", 
"impute","GO.db","preprocessCore","ggplot2", "pheatmap",
"VennDiagram","scatterplot3d","PerformanceAnalytics","org.Hs.eg.db")

chelp=function(){
cat("chelp()")
cat("\nYou can use <install.packages()> to install packages\n")
cat("\nYou can use <library()> to library packages\n")
cat("\n<BiocManager> package is suggested first to install \n")
cat("\n<BiocManager::install> to install some special packages \n")
cat("\nneed packages:\n\n")
cat(paste(cpacks,collapse=";"))
cat("\n")}

cbp=c("#ed1299","#09f9f5","#246b93","#cc8e12","#d561dd","#c93f00","#ddd53e",
"#4aef7b","#e86502","#9ed84e","#39ba30","#6ad157","#8249aa","#99db27","#e07233",
"#ff523f","#ce2523","#f7aa5d","#cebb10","#03827f","#931635","#373bbf","#a1ce4c",
"#ef3bb6","#d66551","#1a918f","#ff66fc","#2927c4","#7149af","#57e559","#8e3af4",
"#f9a270","#22547f","#db5e92","#edd05e","#6f25e8","#0dbc21","#280f7a","#6373ed",
"#5b910f","#7b34c1","#0cf29a","#d80fc1","#dd27ce","#07a301","#167275","#391c82",
"#2baeb5","#925bea","#63ff4f")

gsav=function(x,name,width=8,height=8,dpi=600){
address=getwd() 
if(file.exists("./CDraft")==TRUE){cat("CDraft Existed!\n")}else{
dir.create("./CDraft",recursive=TRUE);cat("CDraft Created\n")}
setwd("./CDraft");ggsave(paste(gsub(":","_",Sys.time()),name,sep="_"),x,width=width,height=height,dpi=dpi)
cat("File had been saved sucessfully under ",paste(address,"/CDraft",sep=""),"\n")
setwd(address)}

gsav2up=function(){
if(file.exists("./CDraft")==TRUE){cat("CDraft Existed!\n")}else{
dir.create("./CDraft",recursive=TRUE);cat("CDraft Created\n")};setwd("./CDraft")}
gsav2down=function(){add=getwd();
cat(paste("File had been saved sucessfully under ",add,"\n",sep=""))}


fsav=function(x,name){
address=getwd() 
if(file.exists("./CWFCSV")==TRUE){cat("CWFCSV Existed!\n")}else{
dir.create("./CWFCSV",recursive=TRUE);cat("CWFCSV Created\n")}
setwd("./CWFCSV")
write.csv(x,paste(name,gsub(":","_",Sys.time()),".csv"),row.names=FALSE)
cat(paste(name,"\nFile had been saved sucessfully under\n",getwd(),"\n",sep=""))
setwd(address)}

pca=function(x){
rownames(x)=x[,1];my_data=x[,-1]
my_data=scale(my_data, center = T, scale = T)
mtcars.pca=prcomp(my_data, center = F,scale. = F)
sm.importance=as.data.frame(summary(mtcars.pca)$importance)
pc1=paste("PC1(",format(sm.importance$PC1[2]*100),"%)",sep="")
pc2=paste("PC2(",format(sm.importance$PC2[2]*100),"%)",sep="")
pc3=paste("PC3(",format(sm.importance$PC3[2]*100),"%)",sep="")
PCA=as.data.frame(mtcars.pca$x)
st=gsub("\\d","",rownames(my_data))
stn=length(levels(factor(st)))
cl=c();for(i in 1:stn){cl=c(cl,rep(cbp[i],length(st)/stn))}
PCA$type =st; PCA$col=cl
return(list(df=PCA,pc1=pc1,pc2=pc2,pc3=pc3))}

cprp=function(x,labxt="Groups",labyt="Proportion%",lname="Type of Germ"){
cann(ncprp);sta=read.csv(deparse(substitute(x)))
prp=ggplot(sta)+geom_bar(aes(x=Group,y=Yv,fill=Bs),stat="identity",position = "fill")+
scale_fill_manual(name=lname,values =sample(cbp,length(levels(factor(sta$Bs))),replace=F))+
theme(panel.background=element_blank(),legend.title=element_text(face ="bold"))+
theme(axis.text.x=element_text(face ="bold"),axis.title=element_text(face ="bold"))+
labs(x=labxt,y=labyt);gsav(prp,"cprp.png");prp}

cbarp=function(x,labxt="Groups",labyt="Values",lname="Types",sp=T,hbar=F){
cann(ncbarp);data=read.csv(deparse(substitute(x)))
if(hbar){pf=data$Yv;shf=c(); shf[pf>=0]="#ff4757"; shf[pf<=0]="#546de5"
barp=ggplot(data,aes(x=Yv,y=reorder(Group,Yv)))+geom_bar(stat="identity",fill=shf)+
theme(axis.title.x = element_text(size = 12,face="bold"),title=element_text(face="bold"),
axis.title.y = element_text(size = 12,face="bold"))+labs(x=labxt,y=labyt)+
theme(panel.background=element_blank(),axis.text.x = element_text(size = 12,color="black"),
axis.text.y = element_text(size = 10,color="black",face="bold"))
gsav(barp,"cbarp.png");return(barp)}
class=levels(factor(data$Group));M=c();S=c()
for (i in 1:length(class)){
M=c(M,mean(data[which(data$Group==class[i]),]$Yv));S=c(S,sd(data[which(data$Group==class[i]),]$Yv))}
barp=ggplot()+geom_bar(aes(x=class,y=M,fill=class),stat="identity",size=1,colour ="black")+
geom_errorbar(aes(x=class,ymin=M-S,ymax=M+S),width=0.4,colour="red",alpha=0.8,size=1)+
labs(x=labxt,y=labyt)+theme(axis.text.x=element_text(face ="bold"),axis.title=element_text(face ="bold"))+
geom_jitter(aes(x=data$Group,y=data$Yv),width=0.1,color="black",size=2)+scale_fill_discrete(name=lname)+
theme(panel.background=element_blank(),legend.title=element_text(face ="bold"))
if(sp){sg=TukeyHSD(aov(Yv~Group,data))$Group; pd=as.numeric(sg[,4][1:length(class)-1])
shsg=c(); shsg[pd>0.05]="no"; shsg[pd<=0.05]="*"; shsg[pd<=0.01]="**"
barp=barp+geom_text(aes(x=class,y=M+S, label =c("",shsg),vjust=0,hjust=0.5),
color="blue4", size=8,show.legend = F,)};gsav(barp,"cbarp.png");barp}

cvlo=function(x,tt=""){
cann(ncvlo);x=read.csv(deparse(substitute(x)))
vlo=ggplot(x)+geom_point(aes(x =LOG,y = -log10(TT), colour=TSS),
alpha=0.4, size=3.5)+scale_color_manual(name="RegTypes",values=c("#546de5", "#d2dae2","#ff4757"))+
geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8)+
geom_hline(yintercept =-log10(0.05),lty=4,col="black",lwd=0.8)+
labs(x="log2(fold change)",y="-log10 (p-value)",title=tt)+
theme_bw()+theme(plot.title=element_text(hjust = 0.5,face ="bold"),legend.position="right",
legend.title=element_text(face ="bold"),axis.title=element_text(face ="bold"))
gsav(vlo,"cvlo.png");vlo}

chot=function(x,w=8,h=8,main="", CR=T,CC=F,SR=T,SC=F,N=F,TYPE=T){
cann(nchot);x=read.csv(deparse(substitute(x)))
c=colorRampPalette(c("green3","black","red3"))(100)
rownames(x)=x[,1];x=x[,-1]
if(TYPE){
getc=colnames(x); st=gsub("\\d","",getc)
annotation_col = data.frame(SampleType =st)
rownames(annotation_col)=getc;lst=levels(factor(st))
if(sum(colnames(x)%in%"Class")>=1){
annotation_row=data.frame(GeneClass=x$Class)
rownames(annotation_row)=rownames(x)
lxc=levels(factor(x$Class));x=x[,-1]
ann_colors = list(SampleType = c(sc=cbp[1:length(lst)]),
GeneClass=c(gr=cbp[(length(lst)+1):(length(lst)+length(lxc))]))
names(ann_colors$SampleType)=lst;names(ann_colors$GeneClass)=lxc
}else{annotation_row=NA;ann_colors = list(SampleType = c(sc=cbp[1:length(lst)]))
names(ann_colors$SampleType)=lst}
cc=seq(length(st)/length(lst), length(st)-length(st)/length(lst),length(st)/length(lst))
hot=pheatmap(x,main=main,scale="row",color=c,cluster_rows=CR, cluster_cols=CC,
border=FALSE,annotation_row=annotation_row,fontface="italic",gaps_col =cc,
fontsize_row=10,fontsize_col =12,annotation_col=annotation_col,show_rownames=SR,
annotation_colors= ann_colors,display_numbers=N,number_color="white")
gsav(hot,"chot.png",w,h)}else{hot=pheatmap(x,main=main,scale="row",color=c,cluster_rows=CR, cluster_cols=CC,
border=FALSE,fontface="italic",fontsize_row=10,fontsize_col =12,show_rownames=SR,show_colnames=SC,
display_numbers=N,number_color="white")};hot}

cven=function(x,t=F,main="",sub=""){
cann(ncven);x=read.csv(deparse(substitute(x)))
vcb=c("red","orange","yellow","blue","violet")
l=list();length(l)=length(colnames(x))
if(length(l)>5){return(cat("Class too much\nThe max cols is five"))}
names(l)=colnames(x)
for(i in 1:length(l)){la=list(unique(x[,i][x[,i]!=""]))
l[i]=la};co=vcb[1:length(l)]
if(t){tt=l;tt$col=co;return(tt)}else{
sname=paste(gsub(":","_",Sys.time()),"cven.PNG",sep="_")
address=getwd();gsav2up();
venn.diagram(l,sname,resolution = 900,alpha=0.5,
fill=co,main=main,main.pos= c(0.5,1.05),cat.cex = 0.7,cat.fontface=4,
main.fontface=4,sub=sub,sub.fontface="bold")
gsav2down();setwd(address)}}

cpca3d=function(x,cex=0.75){
cann(ncpca3d);z=read.csv(deparse(substitute(x)))
PCA=pca(z)$df
address=getwd();gsav2up();
png(paste(gsub(":","_",Sys.time()),"cpca3d.png",
sep="_"),units="in",width=8,height=8,res=600)
with(PCA, scatterplot3d(PC1, PC2, PC3,xlab=pca(z)$pc1,ylab=pca(z)$pc2,zlab=pca(z)$pc3,
pch=16,color=col,cex.symbols = 1.2, font.lab = 2, font.axis = 2))
legend("topleft",legend=levels(factor(PCA$type)),pch=16,col=levels(factor(PCA$col)),ncol=2,cex=cex)
dev.off();gsav2down();setwd(address)}

cpca2d=function(x,tt=""){
cann(ncpca2d);z=read.csv(deparse(substitute(x)))
PCA=pca(z)$df
pca2d=ggplot(PCA)+geom_point(aes(x =PC1,y =PC2, colour=col),
alpha=0.4, size=3.5)+scale_color_discrete(name="Types",
breaks=levels(factor(PCA$col)),labels=levels(factor(PCA$type)))+
labs(x=pca(z)$pc1,y=pca(z)$pc1,title=tt)+
theme_bw()+theme(plot.title=element_text(hjust = 0.5,face ="bold"),legend.position="right",
legend.title=element_text(face ="bold"),axis.title=element_text(face ="bold"))
gsav(pca2d,"cpca2d.png");pca2d}

cwgcna=function(m, n="", main2="Module-trait relationships",
 index3="", main4 = "Network heatmap plot",pw4=2){
cann(ncwgcna);x=read.csv(deparse(substitute(m)))
z=t(x[,-1]);colnames(z)=x[,1];data=as.matrix(z);
colnames(data)=colnames(z);rownames(data)=rownames(z)
powers = c(1:30)
pst=pickSoftThreshold(data, powerVector = powers, verbose = 0)
if(!is.na(pst$powerEstimate)){pw=pst$powerEstimate}else{pw=2}
znet = blockwiseModules(data, power=pw,TOMType = "unsigned", 
minModuleSize=30,reassignThreshold=0, mergeCutHeight=0.25,
numericLabels=TRUE, pamRespectsDendro=FALSE,saveTOMs=TRUE,
saveTOMFileBase="femaleMouseTOM",verbose=3)
mergedColors=labels2colors(znet$colors)
address=getwd();gsav2up();
png(paste(gsub(":","_",Sys.time()),"cWGCNA1.png",
sep="_"),units="in",width=8,height=8,res=600)
plotDendroAndColors(znet$dendrograms[[1]],
cbind(mergedColors[znet$blockGenes[[1]]]),
c("Modules"),dendroLabels=FALSE, hang=0.03,addGuide=TRUE, guideHang=0.05)
dev.off();gsav2down();setwd(address)
moduleLabels=znet$colors; moduleColors=labels2colors(znet$colors)
MEs0=moduleEigengenes(data, moduleColors)$eigengenes
MEs=orderMEs(MEs0)
if(as.character(substitute(n))==""){y=n}else{
y=read.csv(deparse(substitute(n)))}
if(class(y)== "data.frame"){
ysortrn=match(rownames(data), y[,1])
y=y[ysortrn,]; rownames(y)=y[,1];y=y[,-1]
nsamples=nrow(y); moduleTraitCor=cor(MEs,y,use="p")
moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nsamples)
textMatrix=paste(signif(moduleTraitCor, 2),"\n[",signif(moduleTraitPvalue, 1),"]",sep ="")
dim(textMatrix)=dim(moduleTraitCor)
annotation_row=data.frame(Mods=rownames(moduleTraitCor))
rownames(annotation_row)=rownames(moduleTraitCor)
ann_colors = list(Mods=substring(rownames(moduleTraitCor),3))
names(ann_colors$Mods)=rownames(moduleTraitCor)
c=colorRampPalette(c("Blue","white","Red"))(100)
WGCNA2=pheatmap(moduleTraitCor,color=c,cluster_rows=F,cluster_cols=F,
display_numbers=textMatrix,show_rownames=F,border=F,number_color="BROWN",
annotation_row=annotation_row,annotation_colors= ann_colors,
legend_breaks=c(min(moduleTraitCor),min(moduleTraitCor)/2,0,max(moduleTraitCor)/2,max(moduleTraitCor)), 
legend_labels=c("-1.0","-0.5","0","0.5","1.0"))
gsav(WGCNA2,"cWGCNA2.png")
if(index3!=""){MET=MEs;MET$index=y[,index3]
MET=orderMEs(MET);par(mar=c(4,2,1,2),cex=0.9)
address=getwd();gsav2up();
png(paste(gsub(":","_",Sys.time()),"cWGCNA3.png",
sep="_"),units="in",width=8,height=8,res=600)
plotEigengeneNetworks(MET, "",
marDendro=c(0,4,1,2),marHeatmap=c(3,4,1,2),
cex.lab=0.8,xLabelsAngle = 90)
dev.off();gsav2down();setwd(address)}else{
address=getwd();gsav2up();
png(paste(gsub(":","_",Sys.time()),"cWGCNA3a.png",
sep="_"),units="in",width=8,height=8,res=600)
plotEigengeneNetworks(MEs, "",plotDendrograms=F,
marDendro=c(4,4,2,4),cex.lab=0.8,xLabelsAngle = 90)
dev.off();gsav2down();setwd(address)}}
geneTree=znet$dendrograms[[1]]
moduleColors=labels2colors(znet$colors)
dissTOM=1-TOMsimilarityFromExpr(data,power=pw)
plotTOM=dissTOM^pw4; diag(plotTOM)=NA
address=getwd();gsav2up();
png(paste(gsub(":","_",Sys.time()),"cWGCNA4.png",
sep="_"),units="in",width=8,height=8,res=600)
TOMplot(plotTOM,geneTree,moduleColors,main=main4)
dev.off();gsav2down();setwd(address);scol=labels2colors(znet$colors)
names(scol)=names(znet$colors);scot=levels(factor(scol))
for(i in 1:length(scot)){tsco=names(scol[scol==scot[i]])
fsav(tsco,toupper(scot[i]))}
cat("COMPELETED MISSSIONS AND PLESASE CHECK YOUR FILES\n")}

cgsea=function(x,pcf=0.3,sn=10){
cann(ncgesa)
library(ggplot2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
z=read.csv(deparse(substitute(x)))
z=na.omit(z); z=z[!is.infinite(rowSums(z)),]
geneList=z$LOG2;names(geneList)=z$ENTREZID
KEGG_gseresult <- gseKEGG(geneList, pvalueCutoff=pcf)
gsav(dotplot(KEGG_gseresult,showCategory=sn),"cdotplot.png")
for(i in 1:min(dim(KEGG_gseresult)[1],sn)){
gsea=gseaplot2(KEGG_gseresult,i,color="red",base_size = 14, pvalue_table = T)
gsav(gsea,paste(KEGG_gseresult$Description[i],i,".png",sep="_"),height=6)}}

cid=function(f,cfo="",cto="ENTREZID"){
library("org.Hs.eg.db")
det=c("ACCNUM","ALIAS","ENSEMBL","ENSEMBLPROT","ENSEMBLTRANS","ENTREZID","ENZYME","EVIDENCE","EVIDENCEALL" , "GENENAME" ,"GO","GOALL","IPI",  "MAP",  "OMIM", "ONTOLOGY","ONTOLOGYALL", "PATH", "PFAM", "PMID","PROSITE","REFSEQ","SYMBOL","UCSCKG","UNIGENE","UNIPROT" )
if(!as.logical(sum(c(cfo,cto)%in%det))){cat(paste("Your cfo or cto must is one of the follow",det,"\n",sep=""));return()}
x=read.csv(deparse(substitute(f)));x$systemD=1;x[,1]=as.character(x[,1])
y=select(org.Hs.eg.db, keys=x[,1], columns=cto,keytype=cfo)
locay=which(colnames(y)==cto);z=data.frame(NID=y[,locay])
for(i in 1:nrow(x)){z[which(y[,1]%in%x[i,1]),2:ncol(x)]=x[i,-1]}
z=unique(z);colnames(z)=c(cto,colnames(x)[-1]);
z=z[,-which(colnames(x)=="systemD")];fsav(z,"cid_table.csv")}

LOG2=function(x,mx,my){m=c();for( i in 1:nrow(x)){m=c(m,log(mean(unlist(c(x[i,mx])),na.rm=T)/mean(unlist(c(x[i,my])),na.rm=T),2))};m}
PJ=function(x,mx,my){m=c();for( i in 1:nrow(x)){m=c(m,t.test(x[i,mx],x[i,my])$p.value)};m}
WPJ=function(x,mx,my){m=c();for( i in 1:nrow(x)){m=c(m,wilcox.test(as.numeric(x[i,mx]),as.numeric(x[i,my]))$p.value)};m=p.adjust(m,"BH");m}
CTOM=function(c){d=matrix(0,nrow(c),ncol(c));for(i in 1:nrow(c)){for(j in i:ncol(c)){d[i,j]=(sum(c[i,]*c[j,])-3*c[i,j])/(min(sum(c[i,])-1-c[i,j],sum(c[j,])-1-c[i,j])+(1-c[i,j]))}};d}
