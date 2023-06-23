arg <- commandArgs(T)
if(length(arg) < 1){
	cat("Argument: Out_Dir \n")
	quit('no')
}
print(arg)
#arg <- c("/home/wmj/02_Metagenome/16s/JZM_20190119/20190325/Analysis/12_Sig_Distance/Total/Total.KWtest.Sig.difference.txt","/home/wmj/02_Metagenome/16s/JZM_20190119/20190325/bin/Total-grouping.info"   ,"/home/wmj/02_Metagenome/16s/JZM_20190119/20190325/Analysis/19_Sig_PCA/Total"     )


library(vegan)
library(ape) # 用于PCA分析
library(ggplot2) #用于画图
library(ade4);
library(fpc);
outDir <- arg[3];
ff=unlist(strsplit(basename(arg[1]),'\\.txt'))[1]
dis <- read.table(arg[1], header=T,row.names = 1, sep = '\t')
group=read.table(arg[2],header=T,sep='\t', stringsAsFactors = FALSE, comment.char = "")
rownames(group)=as.character(group[,1])
#col.line=c("#DC143C","#4169E1","#2E8B57","#9932CC","#FF8C00","#FFFF00","#8B4513","#FF69B4","#808080")
col.line=unique(as.character(group[,3]))
names(col.line)=unique(as.character(group[,2]))
ii =intersect(rownames(group),rownames(dis))
if(length(ii)!=nrow(group)){
	print('Data Lost')
	ii1=setdiff(rownames(group),ii)
	print(ii1)
}



group=group[ii,]
dis=dis[ii,ii]

Group <- unique(as.character(group[,2])) #按照目的输入样本
shape <-16 #定义点形状
color <- col.line[Group] #定义点颜色
pc=princomp(dis)

Symbol.Group <- as.numeric(as.factor(group[rownames(dis),2]))
Color.Group <- group[rownames(dis),3]
pdf(paste0(outDir,'/',ff,'.PCA.circle.pdf'),width=4,height=4);

plot(pc$loadings[,1],pc$loadings[,2],col = Color.Group,pch=Symbol.Group,xlab=paste("PC1(",sprintf("%.1f",pc$sde[1]^2/sum(pc$sde^2)*100),"%)",spe=""),ylab=paste("PC2(",sprintf("%.1f",pc$sde[2]^2/sum(pc$sde^2)*100),"%)",spe=""),cex=0.9);

s.class(pc$loadings[,1:2], fac=as.factor(group[rownames(dis),2]),grid=F, addaxes=F,axesell =T,label=Group,col=col.line[Group],pch=Symbol.Group,add.plot=T);

dev.off()


pdf(paste0(outDir,'/',ff,'.PCA.no.circle.pdf'),width=4,height=4);

plot(pc$loadings[,1],pc$loadings[,2],col = Color.Group,pch=Symbol.Group,xlab=paste("PC1(",sprintf("%.1f",pc$sde[1]^2/sum(pc$sde^2)*100),"%)",spe=""),ylab=paste("PC2(",sprintf("%.1f",pc$sde[2]^2/sum(pc$sde^2)*100),"%)",spe=""),cex=0.9);

dev.off()

