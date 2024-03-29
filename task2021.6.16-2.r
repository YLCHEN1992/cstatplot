cann=function(x){
up="Map R translation and annotation\n"
down="\nOver 2021.6.16\nneed packages: ggplot2; pheatmap; VennDiagram; scatterplot3d\n"
ncvlo="\ncvlo: need 3 cols \n    LOG->X values (log2);\n    TT->Y values (P-values);\n    TSS->point colors\n"
ncprp="\ncprp: need 3 cols \n    Bs->germ types;\n    Yv->Y values;\n    Group->Treatment type\n"
ncbarp="\ncbarp:  need 2 cols \n    Group->Treatment type; \n    Yv->Y values\n"
nchot="\nchot:  need n cols \n    col1:2->GenenamesID and <Class>; \n    cols...->samplenames\nnote: need rules\n"
ncven="\ncven: you just to fill by cols for items, and the max cols is five. So, It's freedom\n"
ncpca2d="\ncpca2d: need data.frame \n    free.sampleID->sample id need regular;\n    indexn->this is free indexs need dim>2\n"
ncpca3d="\ncpca3d: need data.frame \n    free.sampleID->sample id need regular;\n    indexn->this is free indexs need dim>3\n"
{cat(up);cat(get(deparse(substitute(x))));cat(down)}}

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

chot=function(x,main="", CR=T,CC=F,SR=T,N=F){
cann(nchot);x=read.csv(deparse(substitute(x)))
c=colorRampPalette(c("green3","black","red3"))(100)
rownames(x)=x[,1];x=x[,-1]
getc=colnames(x[,-1]); st=gsub("\\d","",getc)
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
gsav(hot,"chot.png");hot}

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
address=getwd() 
if(file.exists("./CDraft")==TRUE){cat("CDraft Existed!\n")}else{
dir.create("./CDraft",recursive=TRUE);cat("CDraft Created\n")}
setwd("./CDraft")
venn.diagram(l,sname,resolution = 900,alpha=0.5,
fill=co,main=main,main.pos= c(0.5,1.05),cat.cex = 0.7,cat.fontface=4,
main.fontface=4,sub=sub,sub.fontface="bold")
cat("File had been saved sucessfully under ",paste(address,"/CDraft",sep=""),"\n")
setwd(address)}}

cpca3d=function(x,cex=0.75){
cann(ncpca3d);z=read.csv(deparse(substitute(x)))
PCA=pca(z)$df
pca3d=function(){with(PCA, scatterplot3d(PC1, PC2, PC3,xlab=pca(z)$pc1,ylab=pca(z)$pc2,zlab=pca(z)$pc3,
pch=16,color=col,cex.symbols = 1.2, font.lab = 2, font.axis = 2))
legend("topleft",legend=levels(factor(PCA$type)),pch=16,col=levels(factor(PCA$col)),ncol=2,cex=cex)}
cat("Picture need manual Saving");pca3d()}

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



