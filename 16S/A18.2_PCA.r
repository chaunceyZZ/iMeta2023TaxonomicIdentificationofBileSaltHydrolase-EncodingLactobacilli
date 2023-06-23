arg <- commandArgs(T)
if(length(arg) < 1){
	cat("Argument: Out_Dir \n")
	quit('no')
}
print(arg)
#arg <- c("/home/wmj/02_Metagenome/16s/JZM_20190119/20190325/13_betaDiv/Total/unweighted_unifrac_Total.gut.sub2W.filter_4_noChimera_rmAlignFail.otu_table_even12000.pc.txt", "/home/wmj/02_Metagenome/16s/JZM_20190119/20190325/bin/Total-grouping.info" , "/home/wmj/02_Metagenome/16s/JZM_20190119/20190325/Analysis/18_raw_PCA/Total", "2"   )


library(vegan)
library(ape) # 用于PCA分析
library(ggplot2) #用于画图
library(ade4);
library(fpc);


outDir <- arg[3]
grp_col <- arg[4]
ff=unlist(strsplit(basename(arg[1]),'\\.'))[1]

dis <- read.table(arg[1], header=T,row.names = 1, sep = '\t')
groupInfo=read.table(arg[2],header=T,sep='\t', stringsAsFactors = FALSE, comment.char = "")
rownames(groupInfo)=as.character(groupInfo[,1])
#col.line=c("#DC143C","#4169E1","#2E8B57","#9932CC","#FF8C00","#FFFF00","#8B4513","#FF69B4","#808080")
col.line=unique(as.character(groupInfo[,3]))
names(col.line)=unique(as.character(groupInfo[,2]))
ii =intersect(rownames(groupInfo),rownames(dis))

if(length(ii)!=nrow(grpInfo)){
	print('Data Lost')
	ii1=setdiff(rownames(grpInfo),ii)
	print(ii1)
}


groupInfo=groupInfo[ii,]
dis=dis[ii,ii]
pc =princomp(dis)

Symbol.Group <- as.numeric(as.factor(groupInfo[rownames(dis),2]))
Color.Group <- groupInfo[rownames(dis),3]
pdf(paste0(outDir,'/',ff,'.PCA.circle.pdf'),width=4,height=4);

plot(pc$loadings[,1],pc$loadings[,2],col = Color.Group,pch=Symbol.Group,xlab=paste("PC1(",sprintf("%.1f",pc$sde[1]^2/sum(pc$sde^2)*100),"%)",spe=""),ylab=paste("PC2(",sprintf("%.1f",pc$sde[2]^2/sum(pc$sde^2)*100),"%)",spe=""),cex=0.9);

s.class(pc$loadings[,1:2], fac=as.factor(groupInfo[rownames(dis),2]),grid=F, addaxes=F,axesell =T,label=Group,col=col.line[Group],pch=Symbol.Group,add.plot=T);

dev.off()


pdf(paste0(outDir,'/',ff,'.PCA.no.circle.pdf'),width=4,height=4);

plot(pc$loadings[,1],pc$loadings[,2],col = Color.Group,pch=Symbol.Group,xlab=paste("PC1(",sprintf("%.1f",pc$sde[1]^2/sum(pc$sde^2)*100),"%)",spe=""),ylab=paste("PC2(",sprintf("%.1f",pc$sde[2]^2/sum(pc$sde^2)*100),"%)",spe=""),cex=0.9);

dev.off()

