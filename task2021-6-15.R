
cann=function(){
cat("Map R translation and annotation\n")
cat("\ncprp: need 3 cols \n    Bs->germ types; \n    Yv->Y values;\n    Group->Treatment type\n")
cat("cbarp: need 2 cols  \n    Group->Treatment type; \n    Yv->Y values\n")
cat("\nOver 2021.6.15\n")}

cprp=function(x,labxt="Groups",labyt="Proportion%",lname="Type of Germ"){
cann();address=getwd() 
cbp=c("#ed1299","#09f9f5","#246b93","#cc8e12","#d561dd","#c93f00","#ddd53e",
"#4aef7b","#e86502","#9ed84e","#39ba30","#6ad157","#8249aa","#99db27","#e07233",
"#ff523f","#ce2523","#f7aa5d","#cebb10","#03827f","#931635","#373bbf","#a1ce4c",
"#ef3bb6","#d66551","#1a918f","#ff66fc","#2927c4","#7149af","#57e559","#8e3af4",
"#f9a270","#22547f","#db5e92","#edd05e","#6f25e8","#0dbc21","#280f7a","#6373ed",
"#5b910f","#7b34c1","#0cf29a","#d80fc1","#dd27ce","#07a301","#167275","#391c82",
"#2baeb5","#925bea","#63ff4f");sta=read.csv(deparse(substitute(x)))
p=ggplot(sta)+geom_bar(aes(x=Group,y=Yv,fill=Bs),stat="identity",position = "fill")+
scale_fill_manual(name=lname,values =sample(cbp,length(levels(factor(sta$Bs))),replace=F))+
theme(panel.background=element_blank(),legend.title=element_text(face ="bold"))+
theme(axis.text.x=element_text(face ="bold"),axis.title=element_text(face ="bold"))+
labs(x=labxt,y=labyt)
if(file.exists("./CDraft")==TRUE){cat("SSNDraft Existed!\n")}else{
dir.create("./CDraft",recursive=TRUE);cat("SSNDraft Created\n")}
setwd("./CDraft");ggsave(paste(gsub(":","_",Sys.time()),"cprp.png",sep="_"),p,width = 8,height=8,dpi=600)
cat("File had been saved sucessfully under ",paste(address,"/CDraft",sep=""),"\n")
setwd(address);p}

cbarp=function(x,labxt="Groups",labyt="Values",lname="Types"){
cann();address=getwd() 
data=read.csv(deparse(substitute(x)))
class=levels(factor(data$Group))
M=c();S=c()
for (i in 1:length(class)){
M=c(M,mean(data[which(data$Group==class[i]),]$Yv));S=c(S,sd(data[which(data$Group==class[i]),]$Yv))}
p=ggplot()+geom_bar(aes(x=class,y=M,fill=class),stat="identity",size=1,colour ="black")+
geom_errorbar(aes(x=class,ymin=M-S,ymax=M+S),width=0.4,colour="red",alpha=0.8,size=1)+
labs(x=labxt,y=labyt)+theme(axis.text.x=element_text(face ="bold"),axis.title=element_text(face ="bold"))+
geom_jitter(aes(x=data$Group,y=data$Yv),width=0.1,color="black",size=2)+scale_fill_discrete(name=lname)+
theme(panel.background=element_blank(),legend.title=element_text(face ="bold"))
if(file.exists("./CDraft")==TRUE){cat("SSNDraft Existed!\n")}else{
dir.create("./CDraft",recursive=TRUE);cat("SSNDraft Created\n")}
setwd("./CDraft");ggsave(paste(gsub(":","_",Sys.time()),"cbarp.png",sep="_"),p,width = 8,height=8,dpi=600)
cat("File had been saved sucessfully under ",paste(address,"/CDraft",sep=""),"\n")
setwd(address);p}





